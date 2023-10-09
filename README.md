# ninotreve.github.io

本文档列举了POCV在落地过程中遇到的精度问题，包括数值精度问题、算法问题和落地工程问题。注意，其中有少一部分问题尚未完全解决，仍处于 todo 状态。

[TOC]

## 一、float 与 double 变量混合运算所引起的误差
时序计算所有流程用的都是用单精度的 float 变量存储与运算，其中可能与双精度 double 变量打交道的地方在以下三点：

1. 在初始化时序数据时，可能会有以下代码：
```c++
ArrivalTime a = 0.0;
```

以上是赋值时直接赋值了一个 double 对象。也可能直接将 `0.0` 传入函数。如果 `ArrivalTime`（或类似的 `Slew`, `Delay` 等等）是 `float`，这里编译器会做隐式数据类型转换。如果 `ArrivalTime` 是 `std::Variant` 或自定义类，会出现数据类型不匹配，无法赋值的情况。此时应该写作
```c++
ArrivalTime a = 0.0f;
```

2. 在 CalcCellDelay 的各种计算模型中，有可能出现 `double` 字面值，如：
```c++
RetVal CalcCellDelayNLDMRModel::CalcRampWaveform(Task *task) const
{
    auto drVertex = task->GetStage()->GetDrvVertex();
    auto rf = task->rf();
    const auto vh = drVertex->GetSlewUpperThreshold(rf);
    const auto vth = drVertex->GetOutputThreshold(rf);
    const auto vl = drVertex->GetSlewLowerThreshold(rf);
    const auto slewDerateFromLib = drVertex->GetSlewDerateFromLib();
    const auto &rmodelInfo = task->GetRmodelInfo();
    static const Eigen::Array<DefaultScalar, 2, 1> thresholds = {{0.01f, 0.97f}};
    static const Eigen::Array<DefaultScalar, 2, 1> thresholds_shift = thresholds - 0.5f;
    PWLWave<DefaultScalar> pwl_wave{Eigen::ArrayX<DefaultScalar>(rmodelInfo.tf * thresholds_shift + rmodelInfo.td),
                                    Eigen::ArrayX<DefaultScalar>(thresholds)};

    task->SetDriverSlew(rmodelInfo.tf * (vh - vl) / slewDerateFromLib);
    task->SetCellDelay(rmodelInfo.td + rmodelInfo.tf / 2.0 - rmodelInfo.tf * vth);
    task->SetPWLWave(std::move(pwl_wave));
    STA_RETURN_SUCCESS;
}
```
注意 `rmodelInfo.td + rmodelInfo.tf / 2.0 - rmodelInfo.tf * vth` 中，由于 `2.0` 是 `double` 型的字面值，导致 `rmodelInfo.tf / 2.0` 的运算也是用 `double` 计算的。除了 `rmodelInfo.tf * vth` 是 `float` 相乘，整个表达式最后都会用双精度进行运算，并输出 `double` 值。这里，由于 `task->SetCellDelay` 函数仅接受 `float`，发生了隐式类型转换。需要注意精度问题的地方在于，它是先用双精度完成计算后，再通过隐式类型转换变成单精度的。

3. MLDM 查表插值的接口 `table->CoeftableLookupCalc` 用的是双精度。请注意调用时进行的数据类型转换。

## 二、加减法的逆运算

首先，我们来定义加减法逆运算的概念。

在 STA 中，时序数据的计算有先后依赖关系，即先得到各条时序弧上的 Delay，各点 Delay 相加，向前传播得到 ArrivalTime，再将 RequiredTime 向后传播，最后通过 ArrivalTime 和 RequiredTime 的相减得到 Slack。但在代码设计中，因为种种原因需要通过 ArrivalTime 倒推出 Delay，或通过 ArrivalTime 和 Slack 倒推出 RequiredTime，那么这时候就需要用到原来加减法的逆运算。

**定义**：若 $d = \text{InversePlus}(d_1, d_2)$，其中 $d_1$ 的一、二、三阶中心矩分别为 $\mu_1, \sigma_1^2, \gamma_1^3$，$d_2$ 的一、二、三阶中心矩分别为 $\mu_2, \sigma_2^2, \gamma_2^3$，则 $d$ 的一、二、三阶中心矩分别为

$$
\begin{aligned}
\mu &= \mu_1 + \mu_2, \\
\sigma^2 &= |\sigma_1^2 - \sigma_2^2|, \\
\gamma^3 &= \gamma_1^3 + \gamma_2^3.
\end{aligned}
$$

若 $d = \text{InverseMinus}(d_1, d_2)$，其中 $d_1$ 的一、二、三阶中心矩分别为 $\mu_1, \sigma_1^2, \gamma_1^3$，$d_2$ 的一、二、三阶中心矩分别为 $\mu_2, \sigma_2^2, \gamma_2^3$，则 $d$ 的一、二、三阶中心矩分别为

$$
\begin{aligned}
\mu &= \mu_1 - \mu_2, \\
\sigma^2 &= |\sigma_1^2 - \sigma_2^2|, \\
\gamma^3 &= \gamma_1^3 - \gamma_2^3. 
\end{aligned}
$$

例如，从 ${\text{AT}}_{i+1} = {\text{Delay}}_{i} + {\text{AT}}_{i}$，可以得出 ${\text{Delay}}_{i} = \text{InverseMinus}({\text{AT}}_{i+1}, {\text{AT}}_{i})$。又如，从 $\text{Slack} = \text{AT} - \text{RAT}$ (early) 或 $\text{Slack} = \text{RAT} - \text{AT}$ (late)，可以得出 $\text{RAT} = \text{InverseMinus}(\text{AT}, \text{Slack})$ (early) 或 $\text{RAT} = \text{InversePlus}(\text{Slack}, \text{AT})$ (late)。

下面是代码中所有需要用到 `InversePlus` 和 `InverseMinus` 的地方。

1. `src/graph_timer2/timer_calc/path_search/PathEnum.cpp`.

> \[todo\]: 注意这个文件里除了下面的 `PathEnum::MakeBranch` 以外，`PathEnum::MakeUnBranch` 函数也要改，但我们不知道怎么改所以就没改。

```c++
// 根据PathUnit，判断当前vertex和tag对应的所有prev tag，生成新的brch
RetVal PathEnum::MakeBranch(const TimingPath *timingPath, const PathUnit *brchUnit, const PathUnit *prevPathUnit)
{
    ...
    for (auto prevInfo : prevInfoSeq) {
        // prevtagSeq已经包含了当前PathEnd对应的branch，需要去除掉
        ...
        // 计算新的branch（arc分支变化以后）对应的slack值
        Slack lastSlack = timingPath->GetPathEnd()->GetOriginSlack();
        ArrivalTime atInTag = brchUnit->GetArrival();
        Delay arcDelay = TimingPath::GetArcDelay(hdc_, preInfoData.prevEdge, tag, preInfoData.prevArcId);
        Slack brchSlack = 0.0f;
        if (tag->GetMaxMin() == MAX) {
            // 原先是 brchSlack = lastSlack + atInTag - preInfoData.atInPrevTag - arcDelay;
            brchSlack = InversePlus(lastSlack, atInTag) - preInfoData.atInPrevTag - arcDelay;
        } else {
            // 原先是 brchSlack = lastSlack - atInTag + preInfoData.atInPrevTag + arcDelay;
            brchSlack = InverseMinus(lastSlack, atInTag) + preInfoData.atInPrevTag + arcDelay;
        }
        // 生成branch对应的pathend（不包括pathseq）并保存到branch中
        ...
    }
    STA_RETURN_SUCCESS;
}
```

2. `src/graph_timer2/timer_calc/rat/RatTimingCalculator.cpp`:
```c++
bool TimingCalculator::CalcMaxDelayPath(RATInitEnv &ratInitEnv, LaunchArrival *launch, CaptureArrival *capture, TimingChkInfoBase *chkDly)
{
    ...
    ratDataBox_.requiredTime = maxDly + captureAT + subRatClockLatency - ratDataBox_.clkUncertainty - checkDly - CalcPathMargin();
    ratDataBox_.slack = ratDataBox_.requiredTime - launch->LaunchAT(maxDelayPath->IgnoreClkLatency());
    // 原先是 ratDataBox_.propagateRAT = launch->LaunchAT() + ratDataBox_.slack;
    ratDataBox_.propagateRAT = InversePlus(launch->LaunchAT(), ratDataBox_.slack);
    ratDataBox_.dominExcep = exEvalutor_.GetMaxDelayPath();
    checkTimingDone_ = true;
    return true;
}

RetVal TimingCalculator::CalcPathFrom2Clk(RATInitEnv &ratInitEnv, LaunchArrival *launch, CaptureArrival *capture, TimingChkInfoBase *chkDly)
{
    ...
    ratDataBox_.slack = ratDataBox_.requiredTime - pathAT;
    // 原先是 ratDataBox_.propagateRAT = launch->LaunchAT() + ratDataBox_.slack;
    ratDataBox_.propagateRAT = InversePlus(launch->LaunchAT(), ratDataBox_.slack);
    checkTimingDone_ = true;
    STA_RETURN_SUCCESS;
}

RetVal HoldCalculator::CalcPathFrom2Clk(RATInitEnv &ratInitEnv, LaunchArrival *launch, CaptureArrival *capture, TimingChkInfoBase *chkDly)
{
    ...
    ratDataBox_.slack = startAt + launch->LaunchAT() - ratDataBox_.requiredTime;
    // 原先是 ratDataBox_.propagateRAT = launch->LaunchAT() - ratDataBox_.slack;
    ratDataBox_.propagateRAT = InverseMinus(launch->LaunchAT(), ratDataBox_.slack);
    ratDataBox_.startClockOpenEdgeValue = startAt;
    ratDataBox_.endClockCloseEdgeValue = startRat;
    checkTimingDone_ = true;
    ReportCalaculation(ratInitEnv, launch, capture);
    STA_RETURN_SUCCESS;
}
```

3. `src/graph_timer2/timer_task/TimingData.cpp`
```c++
    RetVal TimingData::SaveVertexTimingResults(const GraphTimer2::Vertex *vertex)
    {
        ...
        // (MAX) RAT = AT + Slack
        // 原来是 auto riseRat = riseAT + riseSlack; auto fallRat = fallAT + fallSlack;
        auto riseRat = InversePlus(riseAT, riseSlack);
        auto fallRat = InversePlus(fallAT, fallSlack);
        db::design::TimingData::TimingValueArray rat = {
            std::numeric_limits<float>::quiet_NaN(),
            std::numeric_limits<float>::quiet_NaN(),
            GetBest(riseRat),
            GetBest(fallRat)
        };
        pin->SetRAT(rat, episode_);
        STA_RETURN_SUCCESS;
    }
```

4. `src/graph_timer2/report/report_target/report_cmd/report.cpp`:
```c++
void ReportTiming::InsertBodyDataPathRowInfo(TString& vertexName, unsigned int curIndex, Cfloat launchClockTime) const
{
    ...
    if (curIndex == startPointIndex_) {
        delay = path_[startPointIndex_]->GetVertex()->GetDrvCellDelayByFloat(path_[startPointIndex_]->GetRiseFall(), path_[startPointIndex_]->GetTag()->GetMaxMin());
        csv_delay = delay;
    } else {
        InsertDerateInfo(curIndex, umRowInfo);
        // 非startpoint 判断是否需要插入derate信息   Tfloat GetTimingDerate(MinMax mm, RiseFall rf)
        // 原本是 auto temp_delay = path_[curIndex]->GetArrival() - path_[curIndex - 1]->GetArrival();
        auto temp_delay = InverseMinus(path_[curIndex]->GetArrival(), path_[curIndex - 1]->GetArrival());
        csv_delay = temp_delay;
        // 原本是 delay = path_[curIndex]->GetArrival() - path_[curIndex - 1]->GetArrival();
        delay = GetWorst(InverseMinus(path_[curIndex]->GetArrival(), path_[curIndex - 1]->GetArrival()));
    }
    ...
}

```

## 三、set driving cell 后 loadedDelay - unloadedDelay

见 `src/graph_timer2/timer_calc/delay_calc/CalcMethod.cpp` 中的 `CalcMethod::CalcDrivingCell` 函数。原先的实现是
```c++
Delay drivingCellDelay = loadedDelay - unloadedDelay;
```
Delay 的 mean 应该相减，但 sigma 和 gamma 该以什么方式相减，我们尝试了两种模式，其一是 $\sigma = \sqrt{\sigma_2^2 – \sigma_1^2}$，其二是 $\sigma = \sigma_2 – \sigma_1$, 最后发觉 PT 用的是 sigma 直接相减的方法，因此为了对齐 PT，我们定义了 `DistDirectMinus` 函数，sigma 和 gamma 直接相减。

## 四、平方和立方相减开根导致的数值精度问题

对浮点数进行 $\sqrt{\sigma_1^2 - \sigma_2^2}$ 或 $\sqrt[3]{\gamma_1^3 - \gamma_2^3}$ 操作是一个非常危险的操作。下面一个小实验说明大数相减导致的精度问题：
```c++
    float f1 = 0.23358909;
    float f2 = 0.23358910;
    printf("f1 = %.7f\n", f1);
    printf("f2 = %.7f\n", f2);
    printf("f = %.7f\n", std::sqrt(f2 * f2 - f1 * f1));
```

得到的结果是：
```text
f1 = 0.2335891
f2 = 0.2335891
f = 0.0000610
```
`f1` 和 `f2` 在 float 精度意义下可以认为是相等的，但平方相减后会得到一个小数，这会影响最后的精度。为此，当我们遇到平方差、立方差相减的运算时，我们使用函数 `TruncatedSqrtDiff` 和 `TruncatedCbrtDiff`, 他们的定义如下：
```c++
constexpr float DIST_FLOAT_TRUNCATE_CHECK_TOLERANCE = 1E-5F;

inline float TruncatedDiff(const float a, const float b)
{
    // 如果直接计算a^2 - b^2会放大float的数值误差
    float diff = a - b;
    float maxAB = std::max(std::abs(a), std::abs(b));
    return (DIST_FLOAT_TRUNCATE_CHECK_TOLERANCE * maxAB) < std::abs(diff) ? diff : 0.0f;
}
 
inline float TruncatedSqrtDiff(const float a, const float b)
{
    return std::sqrt(std::abs(TruncatedDiff(a, b) * (a + b)));
}
 
inline float TruncatedCbrtDiff(const float a, const float b)
{
    return std::cbrt(TruncatedDiff(a, b) * (a * a + a * b + b * b));
}
```

## 五、Timing 值的比较
一些场合中可能需要对两个 Timing 值进行比较，常见的有：
+ 通过和 0 比较来确定某个 Timing 为正数。例如 `if (GetDistNominal(inslew) > 0) {...}`.
+ 通过比较两个 Timing 的差值是否小于某个误差界，来判断两个 Timing 是否相同。例如下面在 `LocalEvalChecker::PrintErrorVertex()` 函数的定义内的一段代码，它尝试比较两个 transition 是否相等。

```c++
    for (auto &[vertex, value] : vts_trans_map) {
        if (std::abs(GetDistNominal(value.trans[0]) - GetDistNominal(vertex->GetSlew(RISE, MAX))) > MIN_DIFF) {
            TIMER_LOG_INTERNAL_WARNING("[LocalSanityCheck] error: Graph %d view %s, vertex %s rise TRANS l/g: %f/%f", local_->GetGraphId(), tag.c_str(), vertex->GetPPort().c_str(), GetDistNominal(value.trans[0]), GetDistNominal(vertex->GetSlew(RISE, MAX)));
            result = false;
        }
        if (std::abs(GetDistNominal(value.trans[1]) - GetDistNominal(vertex->GetSlew(FALL, MAX))) > MIN_DIFF) {
            TIMER_LOG_INTERNAL_WARNING("[LocalSanityCheck] error: Graph %d view %s, vertex %s fall TRANS l/g: %f/%f", local_->GetGraphId(), tag.c_str(), vertex->GetPPort().c_str(), GetDistNominal(value.trans[1]), GetDistNominal(vertex->GetSlew(FALL, MAX)));
            result = false;
        }
    }
```

+ 通过比较来选取较为 Worst（或 Best）的值。例如下面一段代码：

```c++
bool PrevTag::MergeFrom(const PrevTag &newPrevTag)
{
    if (newPrevTag.Arc() != nullptr && newPrevTag.Surface() != nullptr) {
        if (this->Similar(newPrevTag)) {
            // 1. timing_type,combinational and combinatinal_fall build 2 path respactively
            // 2. surface type, cell_rise/fall and retaining_rise/fall build 2 path respactively
            // 3. timing_sense, pos/neg + non build 2 path respactively
            if (newPrevTag.IsCheckTimingType() && !IsMergeByTimingType(newPrevTag)) {
                return false;
            }
            if (newPrevTag.IsCheckTimingSense() && !IsMergeByTimingSense(newPrevTag)) {
                return false;
            }
            if (newPrevTag.IsCheckSurfaceType() && !IsMergeBySurfaceType(newPrevTag)) {
                return false;
            }
            if ((this->fromTag_->GetMaxMin() == MAX && GetWorst(newPrevTag.Arrival()) > GetWorst(this->Arrival())) ||
                (this->fromTag_->GetMaxMin() == MIN && GetBest(newPrevTag.Arrival()) < GetBest(this->Arrival()))) {
                this->Update(newPrevTag.Arc(), newPrevTag.Surface(), newPrevTag.Arrival());
                this->fromArcId_ = newPrevTag.fromArcId_;
            }
            return true;
        }
    }
    return false;
}
```

在没有 POCV 的时候，所有的数据都是 float，所有的比较都是自然的。
而当 Timing 变为 variant 类型后，所有的这些比较都应当被明确，到底要比较 Timing 的 `mean`, `mu ± 3sigma` 还是 `nominal`. 明确之后，就应当使用相应的
`GetDistMean`, `GetWorst/GetBest`, `GetDistNominal` 取出相应的值再比较。没有经过这样处理的比较，我们称为 **裸比较**。

裸比较的行为是和 Timing 数值完全无关的。如果直接去比较两个 variant 对象时，C++ 17 会比较两者的类型的 index. 回忆一下，我们定义了
```c++
#define TIMING_VARIANT std::variant<float, GaussDist, NonGaussDist>
```

这个定义决定了 `float`, `GaussDist`, `NonGaussDist` 对应的 index 分别为 0,1,2. 

假设我们进行了如下的裸比较
```c++
TIMING_VARIANT v = 0.0f;
TIMING_VARIANT gv = GaussDist(1.0f, 0.1f);
TIMING_VARIANT ngv = NonGaussDist(1.0f, 1.1f, 0.1f, 0.01f);

bool compare1 = (v > gv);
bool compare2 = (ngv > v);
```
那么 compare1 相当于在比较 `float` 的 index 0 和 `GaussDist` 的 index 1 哪个更大，最终 `compare1 <- False`. compare2 相当于在比较 `NonGaussDist` 的 index 2 和 `float` 的 index 0 哪个更大，最终 `compare2 <- True`. 

因此现在我们完全重载了 `TIMING_VARIANT` 的所有比较运算符，从而使得所有的裸比较都会引发一个 Internal Error. 例如：
```c++
inline bool operator>=(const double&, const TIMING_VARIANT&)
{
    TIMER_LOG_INTERNAL_ERROR("Invalid usage of >= between double and TIMING_VARIANT");
    return false; 
}
```

## 六、关于 GetWorst()/GetMeanAddNSigma() 和 GetBest()/GetMeanMinusNSigma() 的使用场合

> 这是一个仅涉及非高斯计算流程的问题。
>
> \[todo\]：这个问题截止 9.25 我们尚未着手开始解决，仍留有一些写着 `GetWorst()` 的代码应当修改为 `GetMeanAddNSigma()` 和一些 `GetBest()` 应当修改为 `GetMeanMinusNSigma()` ，以后有机会的话我们再来修改这一部分，没机会的话要请华为同事来修改。

首先明确这四个函数的功能如下:
+ `GetWorst(timing_data)`: 返回 `timing_data` 的 99.865% 分位点。对于高斯型随机变量来说，它等价于返回随机变量的 μ+3σ 值. 对于非高斯变量，它通过一个插值模型（通过对 pt 逆向得到）查表得来，不存在简单的等价。
+ `GetMeanAddNSigma(timing_data)`: 返回 `timing_data` 的 μ+3σ 值。
以及
+ `GetBest(timing_data)`: 返回 `timing_data` 的 0.135% 分位点。对于高斯型随机变量来说，它等价于返回随机变量的 μ-3σ 值. 对于非高斯变量，它通过一个插值模型（通过对 pt 逆向得到）查表得来，不存在简单的等价。
+ `GetMeanMinusNSigma(timing_data)`: 返回 `timing_data` 的 μ-3σ 值。

根据我们对 POCV 计算流程的认识，以 `GetWorst()/GetMeanAddNSigma()` 这一对为例，我们认为它们的使用场合的区别是：
+ `GetWorst()`: 用于高斯和非高斯随机变量的 **Report 和 Dump 相关**的场合。例如执行 `rpt_path_timing` 命令时，要打印的 delay 和 at 都应当用 `GetWorst()` 来获取。因为无论是高斯还是非高斯，在 report 其 corner 值的时候，我们直接关心的是 timing 的 99.865% 或 0.135% 分位点。（注意，在 report 或 dump slew 时，应当统一用 `GetDistNominal()` 来转换，详见 **关于 GetDistNominal() 和 GetDistMean() 的使用场合** 一节）. 例如
```c++
// GraphReport.cpp
void GraphReport::DumpAT(const Graph* graph, std::ofstream &out) const
{
    //  vertex,rise_slew,rise_AT,rise_slack,fall_slew,fall_AT,fall_slack,rise_cap,fall_cap, level
    for (const auto &vtx : graph->GetAllVertices()) {
        if (vtx->IsHier()) {
            continue;
        }
        ArrivalReporter maxRiseReporter;
        ArrivalReporter maxFallReporter;
        ArrivalReporter minRiseReporter;
        ArrivalReporter minFallReporter;
        minRiseReporter.mm = MinMax::MIN;
        minFallReporter.mm = MinMax::MIN;
        ArrivalReporter biMaxRiseReporter;
        ArrivalReporter biMaxFallReporter;
        ArrivalReporter biMinRiseReporter;
        ArrivalReporter biMinFallReporter;
        biMinRiseReporter.mm = MinMax::MIN;
        biMinFallReporter.mm = MinMax::MIN;
        Slew slewRise = vtx->GetSlew(RISE, MAX);
        Slew slewFall = vtx->GetSlew(FALL, MAX);
        Slew slewRiseMin = vtx->GetSlew(RISE, MIN);
        Slew slewFallMin = vtx->GetSlew(FALL, MIN);
        auto capRise = vtx->GetCap(RISE, MAX);
        auto capFall = vtx->GetCap(FALL, MAX);
        auto capRiseMin = vtx->GetCap(RISE, MIN);
        auto capFallMin = vtx->GetCap(FALL, MIN);
        GetATGroups(vtx, &maxRiseReporter, &maxFallReporter, &minRiseReporter, &minFallReporter);
        if (vtx->GetDirection() == INOUT) {
            auto biVertex = graph->GetBidirecVertex(vtx);
            if (biVertex == nullptr) {
                continue;
            }
            GetATGroups(biVertex, &biMaxRiseReporter, &biMaxFallReporter, &biMinRiseReporter, &biMinFallReporter);
            maxRiseReporter.MergeFrom(biMaxRiseReporter);
            maxFallReporter.MergeFrom(biMaxFallReporter);
            minRiseReporter.MergeFrom(biMinRiseReporter);
            minFallReporter.MergeFrom(biMinFallReporter);
            capRise = biVertex->GetCap(RISE, MAX);
            capFall = biVertex->GetCap(FALL, MAX);
            capRiseMin = biVertex->GetCap(RISE, MIN);
            capFallMin = biVertex->GetCap(FALL, MIN);
        }
        out << vtx->GetPPort().c_str() << ",";
        out << GetDistNominal(slewRise) << "," << GetWorst(maxRiseReporter.arrivalTime) << "," << GetBest(maxRiseReporter.slack) << ",";
        out << GetDistNominal(slewFall) << "," << GetWorst(maxFallReporter.arrivalTime) << "," << GetBest(maxFallReporter.slack) << ",";
        out << capRise << "," << capFall << ",";
        out << GetDistNominal(slewRiseMin) << "," << GetWorst(minRiseReporter.arrivalTime) << "," << GetBest(minRiseReporter.slack) << ",";
        out << GetDistNominal(slewFallMin) << "," << GetWorst(minFallReporter.arrivalTime) << "," << GetBest(minFallReporter.slack) << ",";
        out << capRiseMin << "," << capFallMin << "," << vtx->GetLevel() << "," << vtx->IsHier()<< std::endl;
    }
}
```

+ `GetMeanAddNSigma()`: 用于高斯和非高斯随机变量的 **Merge 相关**的场合。例如当比较两个 Timing 的大小来确定哪一个更悲观或更乐观。例如
```c++
// GBAData.cpp
bool ArrivalTagData::SaveArrivalTime(MinMax mm, ArrivalTime at, uint8_t source)
{
    // AT Merge using Statistical Graph Pessimism
    ArrivalTime morePessimisticAT;
    bool merged = true;
    if (mm == MAX) {
        morePessimisticAT = MergeDistMax(arrivalTime, at);
        merged = (GetMeanAddNSigma(morePessimisticAT) > GetMeanAddNSigma(arrivalTime));
    } else if (mm == MIN) {
        if (FuzzyEqual(arrivalTime, DEFAULT_ARRIVAL)) {
            arrivalTime = DEFAULT_MIN_ARRIVAL;
        }
        morePessimisticAT = MergeDistMin(arrivalTime, at);
        merged = (GetMeanMinusNSigma(morePessimisticAT) < GetMeanMinusNSigma(arrivalTime));
    } else {
        morePessimisticAT = at;
    }
    if (merged) {
        arrivalBits_.tagSource = source;
    }
    arrivalTime = morePessimisticAT;
    return merged;
}
```

## 七、关于 GetDistNominal() 和 GetDistMean() 的使用场合

> 以下的使用场合的观点没有经过验证，能确定的只有 *primetime 使用 slew 的 nominal 值来查表* 这件事。
>
> \[todo\]: 仍然有大量的 `GetDistMean()` 应该改为 `GetDistNominal()`, `GetWorst()` 和 `GetBest()`，这种情况主要集中在各种 report 和 log 语句中。

顾名思义，`GetDistNominal()` 的功能是获取随机变量的 nominal 值（标称值），`GetDistMean()` 的功能是获取随机变量的 mean 值（平均值）。它们的使用场合为：

1. `GetDistNominal()` 用于所有需要用 delay, slew 的一个 float 值来进行 **与 pocv 没有直接关系的其他计算** 的场合. 例如，计算 Net Delay 时。
2. `GetDistMean()` 仅用于调试和 dump 相关的 variantion 信息，在实际计算流程中不使用。

此外，Slew 的查表用的是 nominal 值而不是 mean 值。Slew Merge 时，不同于 AT Merge 比较 μ±3σ，而 Slew 只比较 nominal 值。以下是现在 slew merge 相关函数的实现。

```c++
bool CalcMethod::SlewMerge(Slew slewTemp, Slew slewNet, MinMax mm)
{
    // return mm == MAX ? FuzzyGreater(slewTemp, slewNet) : FuzzyLess(slewTemp, slewNet);
    return mm == MAX ? FuzzyGreater(GetDistNominal(slewTemp), GetDistNominal(slewNet)) : FuzzyLess(GetDistNominal(slewTemp), GetDistNominal(slewNet));
}

Slew MergeSlew(MinMax mm, Slew s1, Slew s2)
{
    if (mm == MinMax::MAX) {
        return GetDistWithMaxNominal(s1, s2);
    }
    return GetDistWithMinNominal(s1, s2);
}
```

## 八、其他事项

+ `local_eval_checker.cpp` 里有大量的比较两个 timing data 是否相等的代码。在无 pocv 的纯 float 环境下，这些比较是直接的。但在 pocv 下，这些比较到底哪些要用 nominal，哪些
要用 corner 值，还是说应该有其他的比较方式，我们是不确定的。这里需要华为负责这个文件的同事来修改。

+ 目前非高斯是不存储 `mean_shift` 量的，而 `mean_shift` 量又是直接查表得到的。考虑到 `mean_shift` 的量级 (1e-5 或 1e-6) 往往远小于 `mean`，因此从浮点数误差的角度来说，非高斯数据应该存储 (`nominal`, `mean_shift`) 而非存储 (`nominal`, `mean = nominal + mean_shift`)，然后在需要使用 mean 的时候再计算出来。然而，我们当初没有考虑到这个问题，所以实现的是存储 
(`nominal`, `mean`) 的版本。未来如果有需要，华为方面可以修改一下。



## 九、pin slack and path slack
## 十、图悲观

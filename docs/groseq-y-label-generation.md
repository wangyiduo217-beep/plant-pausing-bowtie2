# GRO-seq 链特异性 y 标签构建流程

本页记录植物 pausing 项目从高可信 BAM 构建深度学习 y 标签的实际流程。记录版本为
`groseq_y_article_analog_v2.0_rrna_refiltered`（2026-09-19）。这里的 y 表示**可重复的新生转录区间在固定窗口中的覆盖强度**，不是 pausing index、read count、表达量、peak 概率或单碱基 pause site 真值。

本方法借用了 *Cross-species prediction of histone modifications in plants via deep learning* 的固定窗口建模思路，但输入信号、区间调用器和标签定义均针对 GRO-seq 重新制定。因此应称为 article-analog 标签，不能表述为对该论文标签流程的逐项复现。

## 流程概览

```text
primary MAPQ20 BAM
  -> 按协议保留单端 reads 或双端 R1
  -> 排除已知 rRNA/细胞器坐标
  -> 同一 SRX 的多个 SRR 合并到一个 HOMER tag directory
  -> HOMER findPeaks -style groseq，分别得到正负链转录区间
  -> 每个物种、每条链单独进行 bedtools multiinter
  -> 保留至少 2 个独立 SRX 支持的共识原子区间
  -> 根据 SRX 支持数计算区间置信权重
  -> 1,024 bp 窗口、512 bp 步长滑窗
  -> 分别计算 y_plus 和 y_minus
  -> 排除与 rRNA/细胞器掩膜重叠的窗口
  -> 输出正窗口 y 表
```

## 1. 输入和分析单位

输入为经过下列筛选的坐标排序 BAM：

```bash
samtools view -b -q 20 -F 2820 input.sorted.bam > input.primary_mapq20.bam
samtools index -c input.primary_mapq20.bam
```

`2820 = 4 + 256 + 512 + 2048`，依次排除未比对、secondary、QC-fail 和 supplementary 记录。不按坐标去重，也不排除 duplicate 标志，因为真实 pausing 热点可能产生相同坐标。在进入标签构建前，BAM 和 CSI 均通过存在性检查和 `samtools quickcheck`。

标签构建以 **SRX** 为独立文库单位。一个 SRX 包含多个 SRR 时，所有对应 BAM 一起传入同一个 HOMER tag directory；这些 SRR 只算一个独立支持，避免将技术拆分误计为多个生物学重复。

最终标签纳入 22 个独立 SRX、共 34 个 SRR：

| 物种 | 项目 | 纳入的 SRX（括号内为 SRR） |
|---|---|---|
| 拟南芥 | GSE109974 | SRX3638175 (SRR6661079); SRX3638176 (SRR6661080) |
| 拟南芥 | GSE117014 | SRX4388308 (SRR7518304); SRX4388309 (SRR7518305) |
| 拟南芥 | GSE128698 | SRX5557585 (SRR8767369); SRX5557586 (SRR8767370) |
| 拟南芥 | GSE181488 | SRX11651146 (SRR15347140); SRX11651147 (SRR15347141); SRX11651148 (SRR15347142) |
| 拟南芥 | GSE181598 | SRX11664624 (SRR15362139); SRX11664631 (SRR15362138) |
| 拟南芥 | GSE83108 | SRX1830041 (SRR3647034); SRX1830042 (SRR3647035) |
| 拟南芥 | GSE95301 | SRX2613510 (SRR5313797) |
| 小麦 | GSE178276 v1 | SRX11155086 (SRR14825452–SRR14825457); SRX11155087 (SRR14825458–SRR14825462) |
| 小麦 | GSE178276 v2 | SRX13823254 (SRR17655078–SRR17655079); SRX13823255 (SRR17655080–SRR17655081); SRX13823256 (SRR17655082–SRR17655083) |
| 玉米 | PRJNA788565 | SRX13402909 (SRR17223339); SRX13402910 (SRR17223338); SRX13402911 (SRR17223337) |

拟南芥 SRX2613509 在协议复核中处于 hold 状态，没有进入最终 y 标签。样本筛选限于已确认的野生型、未接受生物学处理的普通 GRO-seq；体外 run-on 条件对照和其他新生转录技术不混入本标签集合。

## 2. 双端 reads 和链方向

GRO-seq 的链信息必须按建库协议解释，不能从 FASTQ 是 R1 还是 R2 直接猜测。

- 单端数据使用通过筛选的全部记录。
- 双端数据只保留 R1，即 `samtools view -f 64`，避免一个片段的两个 mate 被重复计入信号。
- GSE181488 的实证链方向与比对链相反，HOMER 调用时增加 `-rev`。
- 其余纳入文库中，R1 或单端 read 的比对链与 RNA 链一致。
- 正链和负链从区间调用开始一直分开处理，最终分别产生 `y_plus` 和 `y_minus`。

方向规则来自前置协议复核和链特异性 QC。若加入新数据集，必须先重新判断链方向，不能沿用物种默认值。

## 3. rRNA 和细胞器处理

每个物种建立 BED 掩膜，包含已识别的 rRNA 富集坐标，以及参考基因组中可识别的线粒体和叶绿体整条序列。掩膜在两个阶段使用：

1. HOMER 前，从分析 BAM 中排除与掩膜重叠的比对记录；
2. HOMER 后，再从转录区间和最终滑窗中排除与掩膜重叠的项目。

玉米正式版还在 FASTQ 层面先比对 rRNA 参考并过滤，随后重新比对基因组。三个正式玉米输入为 sequence-level rRNA-filtered、primary、MAPQ ≥ 20 BAM。玉米仍沿用相同的坐标掩膜步骤，但本次玉米掩膜没有产生额外排除窗口。

坐标掩膜只能处理已知坐标，序列过滤也受 rRNA 参考完整性限制。两者都不能被描述为对所有 rRNA 污染的绝对清除。

## 4. 每个 SRX 调用链特异性 GRO-seq 区间

同一个 SRX 的一个或多个处理后 BAM 共同建立一个 HOMER tag directory：

```bash
makeTagDirectory <tag_directory> <SRR1.analysis_input.bam> [SRR2.analysis_input.bam ...]
```

随后使用 HOMER 5.1：

```bash
findPeaks <tag_directory> \
  -style groseq \
  -tssFold 4 \
  -bodyFold 3 \
  -minBodySize 500 \
  -maxBodySize 100000 \
  -pseudoCount 1 \
  -o auto
```

仅 GSE181488 在上述命令末尾增加 `-rev`。`transcripts.txt` 经 `pos2bed.pl` 转为 BED，再排除掩膜区间。每个 SRX 的结果都检查 BED 列数、坐标和正负链是否有效。

这里调用的是较宽的新生转录区间，作用类似 interval calling；它不是针对窄富集信号的通用 ChIP-seq peak calling，也不直接定位单碱基 pause site。

## 5. 构建跨文库共识区间

对每个物种的 `+` 和 `-` 链分别执行：

```bash
bedtools sort -i <SRX.strand.bed>
bedtools multiinter -i <SRX1.strand.bed> <SRX2.strand.bed> ... \
  -names <SRX1> <SRX2> ...
```

`multiinter` 将重叠关系切分为不重叠的原子区间，并给出每段由多少个独立 SRX 支持。只保留：

```text
support >= 2 independent SRX
```

同一 SRX 的多个 SRR 已在上一步合并，因此不会重复增加 support。正负链完全独立求交，不把相反链信号合并。

共识 BED 各列为：

| 列 | 含义 |
|---|---|
| 1–3 | `chrom`, `start0`, `end0`，0-based、右端不包含 |
| 4 | 唯一共识区间 ID |
| 5 | `confidence`，范围 0.1–1.0 |
| 6 | `strand`，`+` 或 `-` |
| 7 | `support`，独立 SRX 支持数 |
| 8 | `supporting_srx`，支持该区间的 SRX 名单 |

## 6. 将重复支持数转换为置信权重

对每个物种、每条链单独统计所有保留区间中观察到的最小和最大支持数。区间 \(i\) 的置信权重为：

\[
C_i = 0.1 + 0.9\frac{n_i-n_{\min}}{n_{\max}-n_{\min}}
\]

其中 \(n_i\) 是支持该区间的独立 SRX 数。当该链所有保留区间的支持数完全相同时，设 `confidence = 1.0`。实现中先写成 `support/N` 再在观察范围内缩放；由于同一物种的 \(N\) 是常数，化简后与上式相同。

该权重衡量的是同一物种内部的重复支持强弱。因为各物种的独立文库数不同，0.8 等数值不应直接解释为跨物种校准后的概率。

## 7. 1,024 bp / 512 bp 滑窗计算连续 y

沿每条参考序列从坐标 0 建立固定网格：

- 窗口长度 `W = 1024 bp`；
- 步长 `S = 512 bp`；
- 只保留完整落在参考序列内的窗口；
- 坐标采用 BED 规则 `[start0, end0)`。

对窗口 \(W_j\) 和链 \(s\)，原始分数为：

\[
S_{j,s}=\sum_{i\in s}
\frac{|W_j\cap I_i|}{1024}\times C_i
\]

其中 \(I_i\) 是与窗口重叠的同链共识原子区间，\(C_i\) 是上一节的置信权重。正链和负链分别累计，允许同一窗口同时具有两条链的信号。

原始分数按以下规则写成最终 y：

1. 使用 `floor(score * 10 + 0.5) / 10` 四舍五入到一位小数；
2. 大于 1.0 的值截断为 1.0；
3. 小于 0.1 的值设为 0.0；
4. 若 `y_plus == 0` 且 `y_minus == 0`，该行不写入正窗口表；
5. 与 rRNA 或细胞器掩膜有任意重叠的窗口被排除。

二分类列根据四舍五入前的原始重叠分数生成：

```text
binary_plus  = 1 if raw_plus  > 0 else 0
binary_minus = 1 if raw_minus > 0 else 0
```

因此在另一条链使该窗口得以保留时，理论上可能出现 `binary=1` 但对应 y 因极小重叠被舍入为 `0.0` 的情况。

## 8. y 表字段

输出文件为 `<species>.positive_windows.y.tsv.gz`：

| 字段 | 含义 |
|---|---|
| `species` | 物种名称 |
| `chrom` | 染色体或参考序列名称 |
| `start0`, `end0` | 0-based、右端不包含的 1,024 bp 窗口 |
| `y_plus`, `y_minus` | 正链和负链连续标签，范围 0.0–1.0 |
| `binary_plus`, `binary_minus` | 四舍五入前是否与对应链共识区间发生重叠 |
| `max_support_plus`, `max_support_minus` | 与窗口重叠的对应链区间中的最大独立 SRX 支持数；不是 read 数 |
| `dataset_status` | 标签版本状态 |

拟南芥和小麦使用 `article_analog_primary`；序列层面重新去除 rRNA 后的玉米正式版使用 `article_analog_primary_rrna_refiltered`。

## 9. 当前正式结果

| 物种 | 独立 SRX | SRR | 共识原子区间 | 正窗口 | 掩膜排除窗口 | 状态 |
|---|---:|---:|---:|---:|---:|---|
| *Arabidopsis thaliana* | 14 | 14 | 327,616 | 114,244 | 8 | complete |
| *Triticum aestivum* | 5 | 17 | 359,535 | 336,848 | 147 | complete |
| *Zea mays* | 3 | 3 | 176,201 | 301,479 | 0 | complete, sequence-level rRNA refiltered |
| **合计** | **22** | **34** | **863,352** | **752,571** | **155** | complete |

玉米序列层面去除 rRNA 后重新计算的结果相对临时版减少 1,000 个共识原子区间和 886 个正窗口。新旧版本共有 301,468 个相同坐标的窗口，其中 301,425 个窗口的 `y_plus/y_minus` 完全相同，43 个发生变化。旧结果独有 897 个窗口，新结果新增 11 个窗口。旧玉米结果保留在归档中，正式入口指向 `groseq_y_article_analog_v2.0_rrna_refiltered`。

正式 release 包含摘要 JSON、来源输入身份信息、结果比较、SHA-256 清单和每个物种的共识 BED/y 表；清单校验、gzip 完整性、行数、字段范围及拟南芥/小麦新旧版本字节一致性均已检查。

## 10. 建模前仍需完成的步骤

当前 `positive_windows.y.tsv.gz` 只包含至少一条链具有正标签的窗口，**不包含完整的 y=0 背景集合**。正式训练前还需要：

1. 从不与共识区间、rRNA 和细胞器区域重叠的核基因组窗口中抽取负样本；
2. 明确负样本与正样本的比例、GC 和可比对性匹配规则；
3. 以染色体或同源区块划分训练、验证和测试集，避免相邻重叠窗口及同源序列泄漏；
4. 对重复一致性、文库复杂度、链方向和批次效应做最终 QC；
5. 固定参考基因组版本、窗口序列提取规则和 release 校验值。

模型可以将 `y_plus`、`y_minus` 作为两个连续输出，也可以使用对应 `binary` 列进行双任务二分类。两种目标定义应在训练前固定，并在方法中明确说明。

## 方法边界

- 这些标签描述可重复的新生转录区间，不直接证明启动子近端暂停。
- HOMER `-style groseq` 产生的是转录区间，不能等同于单碱基 pause calling。
- 置信权重来自跨 SRX 支持和窗口覆盖，不是统计显著性或后验概率。
- 1,024 bp 窗口适合当前跨物种序列模型，但会平滑单碱基端点信息。
- 新数据集只有通过协议、链方向和参考版本核对后才能加入现有标签集合。

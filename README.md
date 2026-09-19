# plant-pausing-bowtie2

植物 GRO-seq 的可复用 Bowtie2 处理代码：**SRA/FASTQ → 去接头与协议剪切 → 可选 rRNA 过滤 → 基因组比对 → 排序、CSI 索引与基本质控**。

基于植物 pausing 项目实际运行的分析流程整理。本版提供单样本与顺序批处理入口，支持拟南芥、小麦、玉米等物种的自备参考索引。选择接头和剪切参数的依据是具体建库方法，不能仅凭物种套用配置。生产流程使用 HISAT2 的 GSE181488 应沿用其剪接比对方法。

从高可信 BAM 构建正负链 GRO-seq 训练标签的完整步骤、公式、样本清单、输出字段和当前三物种结果见 [GRO-seq y 标签构建流程](docs/groseq-y-label-generation.md)。该标签流程以 1,024 bp 窗口和 512 bp 步长生成 `y_plus/y_minus`，目前的结果表仅包含正窗口。

## 快速开始

实际分析在 Linux 运行；其他系统可使用 `plan` 检查命令。Python 代码仅使用标准库，生物信息工具通过环境安装。

```bash
conda env create -f environment.yml
conda activate plant-pausing-bowtie2
python -m pip install --no-deps -e .
```

`environment.yml` 固定了生产环境所记录的工具版本。版本来源与验证范围见 [生产来源记录](docs/production-provenance.md)。参考基因组、SRA、FASTQ、BAM 和访问凭据放在各自的数据位置；仓库包含代码、示例配置、文档和测试。

### 1. 准备参考索引

已有完整 Bowtie2 索引时，直接在配置中填写索引前缀。新建索引示例：

```bash
plant-bowtie2 build-index --fasta references/genome.fa --prefix references/genome --threads 8
plant-bowtie2 build-index --fasta references/rrna.fa --prefix references/rrna --threads 8
```

索引前缀不包含 `.1.bt2` 或 `.1.bt2l` 后缀。命令检查六个分片，已有输出时拒绝覆盖。`--dry-run` 可先查看建索引命令。提供 rRNA 索引前，应明确所用 rRNA 序列来源及覆盖范围。

### 2. 按建库协议填写配置

复制一个示例为 `configs/local.json`，修改其中的数据路径、参考前缀、样本名称和协议参数：

| 示例 | 输入与用途 |
|---|---|
| [example_wheat_pe.json](configs/example_wheat_pe.json) | 双端 FASTQ，NEXTflex 配置示例，双端各剪 4 nt 随机碱基 |
| [example_arabidopsis_se.json](configs/example_arabidopsis_se.json) | 单端 FASTQ，接头加 poly(A) 处理示例 |
| [example_sra.json](configs/example_sra.json) | 本地 SRA，声明双端布局后转换并处理 |

示例数据文件名为占位符。所有相对路径均相对于 **JSON 配置所在目录**。每个 `run` 必须唯一；`sample` 写入 BAM 的 `SM` 读组标签，可用实际 SRX 或生物学样本名。

SRA 与 `reads` 两种输入二选一；`reads` 列表含一个单端或两个双端 FASTQ。SRA 输入必须声明 `layout` 为 `SINGLE` 或 `PAIRED`，应先独立完成下载校验。转换发现额外 singleton 时会停止并保留文件，供核对。

### 3. 检查命令，再运行

```bash
plant-bowtie2 plan configs/local.json
plant-bowtie2 run configs/local.json
```

`plan` 展示逐步骤参数，不启动工具、不读取大型输入，也不创建分析输出目录。选择单个样本或若干样本：

```bash
plant-bowtie2 run configs/local.json --sample sample_rep1
plant-bowtie2 run configs/local.json --sample sample_rep1 --sample sample_rep2
```

也可使用 `python -m plant_pausing_bowtie2 ...`。批处理按配置中的样本顺序执行，某个样本失败即停止并报告错误。

## 关键配置

`settings.trimming` 定义共有协议；单个样本也可提供 `trimming` 对象覆盖对应字段。

| 字段 | 含义 |
|---|---|
| `reference_index` | 必填，基因组 Bowtie2 索引前缀 |
| `rrna_index` | 可选，省略时记录为未进行 rRNA 预过滤 |
| `trimming.adapter_r1/adapter_r2` | Cutadapt 3′ 接头列表；表达式作为独立参数传入 |
| `trimming.poly_a` | 是否启用 Cutadapt poly(A) 处理 |
| `trimming.random_clip` | 每端前后各剪的碱基数，默认 0；NEXTflex 示例为 4 |
| `trimming.approved` | 执行前必须为 true；应先核对实际协议 |
| `fastqc` | 默认 true，生成 clean FASTQ 的 FastQC 结果 |
| `compress_raw` | 可选生成经完整性检查的原始 FASTQ gzip 副本，保留输入原件 |

默认线程数：转换 8、Cutadapt 8、rRNA 24、基因组比对 16、排序 4、其他 samtools 4、压缩 8。配置键分别是 `conversion/cutadapt/rrna/alignment/sort/qc/compression`。这些是单样本各阶段的线程数，多个独立进程同时运行时应合计评估 CPU、内存和磁盘需求。

生产服务器采用 4 个在途样本、80 CPU 预算及单转换通道；该资源调度策略的依据见 [生产来源记录](docs/production-provenance.md)。本版入口按样本顺序执行，便于单独复现 Bowtie2 处理步骤。

## 输出与重跑规则

每个样本使用独立目录：

```text
results/<run>/
  complete.json             # 全步骤成功后才写入
  status.json               # 当前阶段或失败原因
  artifacts/
    data/                   # clean/non-rRNA FASTQ，sorted/primary_mapq20 BAM + CSI
    qc/                     # 原始 reads 探针，Cutadapt、flagstat、idxstats、FastQC
    logs/                   # 每步的实际命令和输出
    raw/                    # SRA 转换产物（若使用 SRA 输入）
    raw_archive/            # 可选压缩副本
```

运行期间产物保存在 `.work/`，成功后成为 `artifacts/`。输入、参考、设置与工具版本匹配，且已完成产物未发生可检测改变时，再次运行会跳过该样本。身份检查使用路径、大小和修改时间，不是大文件的全内容哈希校验。部分完成或配置不一致时保留现有文件并报错；核对原因后使用新的输出目录重跑。

样本目录有独占锁；进程失败或中断时会清理其启动的工具进程。此版本保守地拒绝不明覆盖，未提供自动断点续跑。

## 方法边界

基因组比对使用 Bowtie2 `--very-sensitive --end-to-end`；双端增加 `--no-mixed --no-discordant -X 1000 --dovetail`。CSI 支持小麦的长染色体。高可信 BAM 使用 `MAPQ >= 20` 和排除标志 `2820`，保留重复标记与细胞器比对记录。

本流程不做坐标去重。rRNA 的双端预过滤仅移除 concordant rRNA pairs，并非穷尽的污染去除。比对和基本质控完成后，仍需检查链方向、复杂度、生物学重复一致性与 pausing 标签，才能评估是否适合作为模型输入。详细参数及项目特异规则见 [分析方法](docs/methods.md)。

当前标签构建方法及建模前尚需补充的负样本步骤见 [GRO-seq y 标签构建流程](docs/groseq-y-label-generation.md)。

## 测试

```bash
python -m unittest discover -s tests -v
python tests/smoke_real_tools.py
```

第一条运行配置与命令构建等测试；第二条需要已安装的真实生物信息工具，使用隔离的小型合成参考和 reads。小样本测试用于验证代码行为，不能替代真实研究样本的全量质控。

本版已通过 Linux 17 项测试和真实双端、单端小数据流程检查；具体结果及未覆盖范围见 [软件验证记录](docs/testing.md)。

## 工具参考

- [Bowtie2 官方手册](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
- [NCBI fasterq-dump 说明](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump)
- [Cutadapt 官方文档](https://cutadapt.readthedocs.io/en/stable/)
- [Samtools 官方文档](https://www.htslib.org/doc/samtools.html)

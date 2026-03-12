# RNA Assessment

RNA 三维结构评估代码库，整理后提供四类稳定能力：

- PDB 归一化
- 基于索引片段的 RMSD 与 P-value 计算
- 基于 MC-Annotate 输出的 Interaction Network Fidelity / Deformation Index
- 纯 Python 的 all-atom lDDT 计算
- 自动链/残基映射与批量 benchmark

这次重构主要解决了原仓库的几个核心问题：源码、示例、构建产物和第三方二进制混放；Python 2/3 语法混杂；外部工具在导入阶段就被强绑定；README 声明的部分能力并没有形成可运行的软件接口。

## 安装

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
```

安装完成后可以直接用 CLI：

```bash
rna-assessment assess \
  examples/data/14_solution_0.pdb \
  examples/data/14_ChenPostExp_2.pdb \
  --per-residue
```

如果同目录存在 `*.index`，工具会自动发现并使用；如果没有，则会根据链级序列对齐自动推断可比较的残基映射。

如果只想算某一个指标，也可以分别使用：

```bash
rna-assessment inf \
  examples/data/14_solution_0.pdb \
  examples/data/14_solution_0.index \
  examples/data/14_ChenPostExp_2.pdb \
  examples/data/14_ChenPostExp_2.index

rna-assessment lddt \
  examples/data/14_solution_0.pdb \
  examples/data/14_ChenPostExp_2.pdb \
  --native-index examples/data/14_solution_0.index \
  --prediction-index examples/data/14_ChenPostExp_2.index

rna-assessment map \
  examples/data/14_solution_0.pdb \
  examples/data/14_ChenPostExp_2.pdb

rna-assessment benchmark \
  examples/data/14_solution_0.pdb \
  examples/data/14_solution_0.pdb \
  examples/data/14_ChenPostExp_2.pdb \
  --native-index examples/data/14_solution_0.index \
  --sort-by rmsd

rna-assessment benchmark \
  --manifest benchmark_manifest.json \
  --sort-by lddt \
  --descending \
  --per-residue
```

`benchmark_manifest.json` 或 `benchmark_manifest.csv` 支持这些字段：

- `prediction` 或 `model`：必填
- `native` 或 `reference`：可选，未提供时可由命令行 `benchmark <native>` 统一补
- `label`
- `native_index`
- `prediction_index`
- `native_annotation`
- `prediction_annotation`

JSON 示例：

```json
[
  {
    "label": "method_a",
    "native": "examples/data/14_solution_0.pdb",
    "native_index": "examples/data/14_solution_0.index",
    "prediction": "examples/data/14_ChenPostExp_2.pdb"
  }
]
```

CSV 示例：

```csv
label,native,native_index,prediction,prediction_index
method_a,examples/data/14_solution_0.pdb,examples/data/14_solution_0.index,examples/data/14_ChenPostExp_2.pdb,
```

也可以用 Python API：

```python
from pathlib import Path

from rna_assessment import calculate_assessment, normalize_structure, prepare_structure_pair

data_dir = Path("examples/data")
output_dir = Path("examples/output")
output_dir.mkdir(exist_ok=True)

normalize_structure(
    data_dir / "14_solution_0.pdb",
    output_dir / "14_solution_0.normalized.pdb",
)

assessment = calculate_assessment(
    data_dir / "14_solution_0.pdb",
    None,
    data_dir / "14_ChenPostExp_2.pdb",
    None,
    include_per_residue=True,
)
print(assessment)

prepared = prepare_structure_pair(
    data_dir / "14_solution_0.pdb",
    None,
    data_dir / "14_ChenPostExp_2.pdb",
    None,
)
print(prepared.native_index, prepared.prediction_index)
```

如果你的 `.mcout` 文件不在 PDB 同目录，可以显式传入：

```bash
rna-assessment assess ref.pdb pred.pdb \
  --native-index ref.index \
  --prediction-index pred.index \
  --native-annotation ref_custom.mcout \
  --prediction-annotation pred_custom.mcout
```

## 项目布局

```text
src/rna_assessment/    核心包
src/RNA_normalizer/    旧包名兼容层
examples/data/         样例输入与预计算 .mcout
third_party/bin/       可执行二进制
third_party/lib/       可选 jar 包
tests/                 自动化测试
```

## 第三方工具策略

核心流程只依赖 `biopython`。

`MC-Annotate` 只在计算 INF 且缺少同目录 `.mcout` 文件时才会被调用。仓库保留了历史 Linux 二进制在 `third_party/bin/MC-Annotate`，这对 Linux 环境可直接用；在 macOS 上通常不可执行，因此示例和测试默认复用仓库内预计算的 `.mcout` 文件，不再要求用户先手工改源码路径。

如果你只给 `PDB` 文件而没有给 `.index`，工具会优先查找同目录 sidecar index，例如 `model.pdb -> model.index`；找不到时再自动推断链映射和可比较残基。

`lDDT` 现在由仓库内部直接计算，不再依赖 OpenStructure。当前实现基于匹配原子上的参考邻域距离，默认使用 15A inclusion radius 与 0.5/1/2/4A 阈值。

如果传入 `--per-residue` 或 `include_per_residue=True`，结果里还会包含每个匹配残基的：

- `lddt`
- `local_rmsd`
- `mean_absolute_error`
- `max_absolute_error`
- `matched_atoms`

`MCQ` jar 被保留在 `third_party/lib/`，但运行它仍然需要本机有 Java Runtime。`GDT`、`ARES`、`MolProbity` 在旧仓库里仍未形成可维护的集成实现，因此当前 CLI 的稳定范围集中在 RMSD、INF 和 lDDT。

## 示例

可以直接运行：

```bash
python examples/basic_usage.py
```

也可以换成自己的输入：

```bash
python examples/basic_usage.py ref.pdb pred.pdb \
  --reference-index ref.index \
  --prediction-index pred.index \
  --per-residue
```

## 测试

```bash
pytest
```

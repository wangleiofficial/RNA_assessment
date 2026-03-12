# RNA Assessment

RNA 三维结构评估代码库，整理后提供三类稳定能力：

- PDB 归一化
- 基于索引片段的 RMSD 与 P-value 计算
- 基于 MC-Annotate 输出的 Interaction Network Fidelity / Deformation Index

这次重构主要解决了原仓库的几个核心问题：源码、示例、构建产物和第三方二进制混放；Python 2/3 语法混杂；外部工具在导入阶段就被强绑定；README 声明的部分能力并没有形成可运行的软件接口。

## 安装

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
```

安装完成后可以直接用 CLI：

```bash
rna-assessment rmsd \
  examples/data/14_solution_0.pdb \
  examples/data/14_solution_0.index \
  examples/data/14_ChenPostExp_2.pdb \
  examples/data/14_ChenPostExp_2.index
```

也可以用 Python API：

```python
from pathlib import Path

from rna_assessment import calculate_interaction_network_fidelity, calculate_rmsd, normalize_structure

data_dir = Path("examples/data")
output_dir = Path("examples/output")
output_dir.mkdir(exist_ok=True)

normalize_structure(
    data_dir / "14_solution_0.pdb",
    output_dir / "14_solution_0.normalized.pdb",
)

rmsd = calculate_rmsd(
    data_dir / "14_solution_0.pdb",
    data_dir / "14_solution_0.index",
    data_dir / "14_ChenPostExp_2.pdb",
    data_dir / "14_ChenPostExp_2.index",
)
print(rmsd)

inf = calculate_interaction_network_fidelity(
    data_dir / "14_solution_0.pdb",
    data_dir / "14_solution_0.index",
    data_dir / "14_ChenPostExp_2.pdb",
    data_dir / "14_ChenPostExp_2.index",
)
print(inf)
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

`MCQ` jar 被保留在 `third_party/lib/`，但运行它仍然需要本机有 Java Runtime。`GDT`、`lDDT`、`ARES`、`MolProbity` 在旧仓库里并没有形成可维护的集成实现，这次重构没有继续保留伪入口，而是把范围收敛到当前真正可验证的核心指标。

## 示例

可以直接运行：

```bash
python examples/basic_usage.py
```

## 测试

```bash
pytest
```

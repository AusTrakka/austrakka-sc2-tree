# austrakka-sc2-tree

[![PyPI - Version](https://img.shields.io/pypi/v/austrakka-sc2-tree.svg)](https://pypi.org/project/austrakka-sc2-tree)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/austrakka-sc2-tree.svg)](https://pypi.org/project/austrakka-sc2-tree)

-----

**Table of Contents**
- [Documentation](https://austrakka.github.io/austrakka-sc2-tree/)
- [Installation](#installation)
- [CLI](#cli)
- [Pipeline](#pipeline)
- [License](#license)

## Installation

```console
pip install austrakka-sc2-tree
```

## Usage

```bash
austrakka-sc2-tree run \
    --fasta tests/data/test.fasta \
    --data tests/data/test.metadata.csv \
    --days-ago 200 \
    --outdir results \
    --name sc2 \
    --tree-threads 8 \
    --dated
```

## CLI 

![](docs/images/cli.png)

## Pipeline

![](docs/images/dag.png)

## License

`austrakka-sc2-tree` is distributed under the terms of the [MIT](https://spdx.org/licenses/MIT.html) license.

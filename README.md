# CDK: b.next Cell Developer Kit
*Internal use only*

The CDK is a set of tools and libraries for building synthetic cells, the operations stack for synthetic cell engineering. This package is the core library, modular code that can be used in specific experiments or in other applications.

## Features
Currently, the CDK contains our core analysis functionality: plate reader and liposome analysis.

## Guides
- [Discovery plate concentration pipeline](docs/discovery_plate_pipeline.md)

## Installation
*TBD*

## Development
### Install poetry
The CDK uses poetry for dependency control and packaging. Install poetry, and activate it to download the dependencies. You can use poetry to manage the development virtual environment (recommended), or create a new conda environment to develop in.

#### Linux
To install poetry on Linux, you can use the following command:

```bash
curl -sSL https://install.python-poetry.org | python3 -
```

#### Mac
*(untested)* Install poetry using homebrew:
```bash
brew install poetry
```

### Activate poetry and download dependencies
```bash
poetry install
```

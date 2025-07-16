# CDK: Nucleus Cell Developer Kit
The Nucleus Cell Developer Kit (CDK) is a set of tools and libraries for building synthetic cells, the stack for synthetic cell engineering. This package is the core library, modular code that can be used in specific experiments or in other applications.

Nucleus is an open-source project for building and working with synthetic cells. For more information, please check out the [Nucleus documentation](https://nucleus.bnext.bio/).

## Features
Currently, the CDK contains our core analysis functionality: plate reader and liposome analysis. In particular, the plate reader code is designed to allow for easy analysis of kinetic timeseries of PURE experiments from Agilent/Biotek and Revvity Envision plate readers. The liposome analysis code is under heavy development, and the code included here is somewhat out of date---please contact us for more information.

## Installation
`pip install nucleus-cdk`

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

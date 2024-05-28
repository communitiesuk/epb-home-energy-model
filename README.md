![test status](https://github.com/communitiesuk/epb-home-energy-model/actions/workflows/test.yml/badge.svg)

# Home Energy Model (HEM) engine written in Rust

## Overview

The [Home Energy Model (HEM)](https://www.gov.uk/government/publications/home-energy-model-technical-documentation) is a
methodology for calculating the energy performance of buildings that is currently under development by BRE under
instruction from the UK Department for Energy Security and Net Zero (DESNZ).

This project constitutes a port into Rust of
the [calculation engine](https://dev.azure.com/BreGroup/_git/Home%20Energy%20Model) written in Python produced by BRE as
a candidate specification of the Home Energy Model.

Its purpose is to provide a library that is significantly more performant than the Python specification code, as well as
providing other properties to support its running as part of a live service, including well-defined errors and
instrumentation for observability.

The project is as yet functionally incomplete, but we have begun working on a test harness that will (in time) be able
to demonstrate complete 1:1 behavioural parity with the specification code.

## Running the engine

Requires the `rustup` toolchain to use. ([Instructions](https://rustup.rs) for installation. For macOS, don't use
Homebrew to install Rust.)

To run tests:

```bash
cargo test
```

Fully running the engine requires a weather file in EPW format and an input JSON file (format not yet documented/
stable):

```
cargo run --release --features="clap indicatif" -- path/to/input.json -e path/to/weather_file.epw
```

The `clap` feature above is mandatory. The `indicatif` feature switches on the output of a progress bar.

### Reverse-engineering JSON schema from input definition

There is currently a package in the project that is able to reverse-engineer a(n incomplete) JSON Schema document from
the HEM input as currently mapped in the Rust code.

Running the following will output this JSON Schema document to stdout:

```bash
cargo run -p schema-gen --release
```

It is expected that a canonical JSON Schema document representing inputs (and outputs) will be published in the near
future, at which point this package will become redundant,
and will be removed.

## Contributing

### Using the commit template

If you've done work in a pair or ensemble why not add your co-author(s) to the commit? This way everyone involved is
given credit and people know who they can approach for questions about specific commits. To make this easy there is a
commit template with a list of regular contributors to this code base. You will find it at the root of this
project: `commit_template.txt`. Each row represents a possible co-author, however everyone is commented out by default (
using `#`), and any row that is commented out will not show up in the commit.

#### Editing the template

If your name is not in the `commit_template.txt` yet, edit the file and add a new row with your details, following the
format `#Co-Authored-By: Name <email>`, e.g. `#Co-Authored-By: Maja <maja@gmail.com>`. The email must match the email
you use for your GitHub account. To protect your privacy, you can activate and use your noreply GitHub addresses (find
it in GitHub under Settings > Emails > Keep my email addresses private).

#### Getting set up

To apply the commit template navigate to the root of the project in a terminal and
use: `git config commit.template commit_template.txt`. This will edit your local git config for the project and apply
the template to every future commit.

#### Using the template (committing with co-authors)

When creating a new commit, edit your commit (e.g. using vim, or a code editor) and delete the `#` in front of any
co-author(s) you want to credit. This means that it's probably easier and quicker to use `git commit` (instead
of `git commit -m ""` followed by a `git commit --amend`), as it will show you the commit template content for you to
edit.

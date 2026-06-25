# Welcome to the CAM Development Repository

This is the repository for CAM documentation. The latest released version of the docs is [v6.3](https://ncar.github.io/cam/cam6.3-users-guide/index.html).

For CAM release code and public (AMWG) development, see the [ESCOMP/CAM](https://github.com/ESCOMP/CAM) repository.

# Modifying the CAM documentation
## Setup
1. If you have not already done so, create a fork of this repository (note that you will likely need to rename the fork because CAM will probably be taken by your fork of ESCOMP/CAM)
1. Clone your forked repository.
1. Either sync your fork to the latest from NCAR (you'll want the latest for the guide/branch you're updating!), or add NCAR as a remote.
1. Branch off of the guide you are updating. The source branch will most likely be `main` (if you're updating the latest development guide), but will be `cam6-users-guide` for CAM6, for example.

If you've synced your fork and wish to create a branch off of main (latest development) called `doc-updates`, the commands are below:
```
git branch doc-updates main
git checkout doc-updates
git submodule update --init (initializes the doc-builder submodule)
```

## Make & commit changes
1. Modify the files you wish to modify (all source files live in `doc/source/users_guide`)
2. Commit the changes to your fork.
```
git add <file changed>
git commit -m '<commit message>'
git push origin doc-updates
```

## Confirm that the docs build
There are several options for building and viewing the docs that you have modified. You can find those options [here](https://escomp.github.io/CTSM/users_guide/working-with-documentation/index.html), though building on casper is recommended.

## Make a PR
Open a PR to the branch/version you are updating.

# Additional notes and tools
## Reviewing a PR
To ease the review process, it is recommended to view the rendered docs rather than reviewing the raw .rst and .md files. To do so, run the following commands in a clone of this repo or or your fork:
```
cd doc/doc-builder/tools
./preview-docs-pr <PR URL>
```

The script will clone the PR repo in your scratch directory (or you can specify a directory elsewhere with the `--dir` flag). Additionally, the script will log information about what files were modified by the PR!

Follow the instructions from the [documentation documentation](https://escomp.github.io/CTSM/users_guide/working-with-documentation/index.html) to view the built docs.

## Bringing in updated doc-builder
When a new tag is made in the [doc-builder repository](https://github.com/ESMCI/doc-builder), follow the steps above to make changes to the documentation (checking out a branch of your fork, etc), and then:

```
git submodule update --remote
git add doc/doc-builder
git commit -m 'update doc-builder submodule'
git push origin <branch name>
```

Then create a PR to the branch in question.

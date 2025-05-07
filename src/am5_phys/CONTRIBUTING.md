# Contributing To am5_phys

Thank you for taking time to contribute.

am5_phys is the repository where AM5 phys code exists. 

What follows is a set of guidelines and how-tos for contributing to am5_phys.
These are guidelines, not rules.  Use your best judgement and feel free to
propose changes to this document in a merge request.

Table of Contents
* [Quick Start Workflow](#quick-start-workflow)
* [Policies](#policies)
* [Merge Requests](#merge-requests)
* [Reporting bugs or issues](#Reporting bugs or issues)
* [Tests](#tests)
* [Style Guide](#styleguides)
* [Release Schedule](#release-schedule)

## Quick Start Workflow

To contribute to this repository see this [step by step guide](MERGE_FOR_AM5_PHYS.md)

## Policies 

* Branch names should be descriptive of what they change.  A branch name of `bugFix` is not specific enough
* Commit messages must be descriptive.  `"Typo"` or `"Bug fix"` is not an acceptable git message. Write the
message with the knowlege that someone will look at it in a year and want to know what this commit does.
* Commit messages should be in the present tense. `"Updates the coupler_nml hyperthreading switch"` 
* Branches must be kept up-to-date with the main branch.

## Merge Requests

Submit merge requests for bug fixes, improvements, including tests, or alerts to
particular problems in the code.  We perform merge requests based on an internal
GFDL schedule that addresses the needs of the GFDL scientists.  This release
schedule is better described in the Release Schedule sections.  

Please keep the changes in a single merge to be as small as possible to help
reviewer(s) quickly evaluate changes.  If you do have a large update, try to
split the update into small, logical merge requests.

Once a merge request is created, a maintainer of the am5_phys repository will 
review the changes, and, if necessary, will work with the author of the merge 
request to modify their code changes. Note that merging merge requests is
contingent on several factors, including the ability to support the changes
long-term, portability, and the scope of the impact on the code base. Therefore,
we do not guarantee that all merge requests will be accepted,
even if the changes pass the initial testing phases, and are otherwise correct.

## Reporting bugs or issues
[Open a new issue using the web interface](https://gitlab.gfdl.noaa.gov/fms/am5_phys/-/issues/new) describing the bug you are solving or the feature you are adding to the code. 

- The issue title should be short and descriptive. 
- The issue description should be clear and concise. Include enough information to help others reproduce the issue, or understand the change requested. 
- Assign the issue to the person that will be fixing it.

## Tests

Users are reponsible for testing their updates and ensuring reproducibility.  If 
reproducibility is broken, then the user must document the reasons for changing the
answers, and this must be accepted by a maintainer of the am5_phys.

If adding new experiments, regression tests and continuous integration must be updated
 to support the new experiment.  More on this as we work through what this will look like.

## Style Guide
Code updates should follow the coding style for the project, contained in the [style guide](#STYLE.md)

## Release Schedule

Releases will be tagged using the format yyyy.rr[.pp], where yyyy is the 4-digit year, rr is the 2-digit release number, and pp is the 2-digit patch number. Preliminary releases mean for testing (i.e., code that is still under development) will be marked yyyy.rr.alpha.pp or yyyy.rr.beta.pp. Alpha tags mark code updates that are intended for developers to include in their baseline regression tests to determine whether the code contains bugs not encountered during baseline testing. Beta tags are intended for a wider audience of developers and end users to ensure that their simulations run as expected and reproduce prior results.

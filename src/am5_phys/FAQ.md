# Git FAQ and Best Practices

The gitlab git cheat sheet is available [here](https://about.gitlab.com/images/press/git-cheat-sheet.pdf).
It can be a useful reference for commands, as well as including install directions and diagrams.
These instructions assume you have already created a fork of this repo via the gitlab site.
If you have not created a fork yet, instructions are included as the first steps in the [merge instructions](MERGE_FOR_AM5_PHYS.md).

### Remotes and Forking

**How do I set up a local repository for development?**
* For general git usage, it is best to clone and have remotes set up for both your fork (https://gitlab.gfdl.noaa.gov/gitlab_username/am5_phys) and
the upstream repo (https://gitlab.gfdl.noaa.gov/fms/am5_phys). This allows you to pull in updates from the upstream repo in order for your code to be up to date with any recent changes.
If you have not created a fork yet, instructions are included as the first steps in the [merge instructions](MERGE_FOR_AM5_PHYS.md).
By default you will have one remote named `origin` pointed to wherever you originally cloned from, which can be checked with `git remote -v`.

After you run `git clone`, to remove the default remote 'origin':
```
git remote rm origin
```

Then add remotes for your fork and the upstream repo:
```
git remote add upstream https://gitlab.gfdl.noaa.gov/fms/am5_phys.git
git remote add myfork https://gitlab.gfdl.noaa.gov/{gitlab_username}/am5_phys.git
```

**How to regularly update your main branch in your forked repository?**
* When you create a fork it will copy the `main` branch from the fms/am5_phys gitlab to your fork, but
after that any changes to the `main` in fms/am5_phys will not go into your forks `main` branch.
You will need to manually pull the changes and then push them in order for it to be up to date.
It is considered a best practice to keep your fork's main up to date, and should always be done
before starting a new branch in order to mitigate potential merge conflicts (conflicting changes to the same file locations).

To update your forks's main from a fresh clone:
```
git clone https://gitlab.gfdl.noaa.gov/<gitlab username>/am5_phys.git
cd am5_phys
git remote add upstream https://gitlab.gfdl.noaa.gov/fms/am5_phys.git
git pull upstream main
git push origin main
```
You will need to enter your gitlab user name and password for the pull and push commands, but after
the push you should see the up-to-date main on your fork's web page.

If want to update your fork's main from an existing local repository, make sure it has remotes
for your fork (<gitlab username>/am5_phys) and the upstream repository (fms/am5_phys). The remote names
may differ depending on how the repo was configured, so you can use the `git remote -v` to ensure
you pull from the fms/am5_phys url and push to your forked copy.
Then you can run the same pull and push commands:
```
git remote -v     # should show both fms/am5_phys url and your fork's url, otherwise see above
git checkout main # switch to current local main branch
git pull upstream main
git push myfork main
```

**Where to push new branches?**
* Any new branches created should be pushed to your fork ({gitlab_username}/am5_phys). They should
also be created off of an up-to-date main branch (see above).

**How do I pull in any new updates to development branches?**
* Once you create a branch, updates may come after the time of creation that will need to be brought
in before it is merged. To bring these updates in, they will need to be fetched (downloaded) and
then 'merged' into your branch to combine the changes your code.

You should push your any changes on your branch to your fork before following the steps below to be safe.

Sometimes the changes can be merged in easily, but if your branch changes code that had also bee
changed by someone else it will lead to merge conflicts that will need to be manually resolved.
For how to fix conflicts, see [Resolving merge conflicts](MERGE_FOR_AM5_PHYS.md#resolving-merge-conflicts).

To fetch main and merge it into a branch:
```bash
## after using git add, commit and push to add changes to YOURbranch
git remote add upstream https://gitlab.gfdl.noaa.gov/fms/am5_phys.git # if not already present
git fetch upstream main
git checkout YOURbranch
# Resolve any conflicts
git merge upstream/main
git push origin YOURbranch
```

**How to share my changes with someone else?**
* When sharing code updates, users should share a branch on their fork and avoid sending a local
path to another user.  By sharing using git, we can ensure that everyone's code is tracked on git and
the history is traceable, and everyone follows a standard for sharing.
Make sure you add, commit (with a descriptive commit message) and push your updates to your fork's remote.
Assuming `origin` is your forked repository:
```bash
git checkout -b changesToShare #make your branch name descriptive
git add -u
git commit "Present tense descriptive message about your updates"
git push origin changesToShare
```
You can see these updates by adding a remote to the person sharing's fork, then fetching and checking out the other persons branch if you are the one
receiving the updates.
```bash
git remote add theirFork https://gitlab.gfdl.noaa.gov/{their_gitlab_username}/am5_phys.git
git fetch theirFork changesToShare
git checkout theirFork/myBranch
```

#### Branches and Local Changes

**How to checkout a new branch?**
* The `main` branch of this repository is protected and can only be updated by people with at least *Maintainer*
privilege.  Users can create user branches in their forked repositories and make their code changes. Use a descriptive name for your branch.
Branch names like `bugFix` or `myBranch` or `Updates` are not useful.
```bash
git checkout -b c192RadNmlUpdate
```

**How to submit a merge request?**
* If a user would like to update the main branch in this repository it can be done
by submitting a Merge Request on the gitlab site https://gitlab.gfdl.noaa.gov/fms/am5_phys/-/merge_requests/new .
More information on submitting merge requests can be found in the [merge request documentation](MERGE_FOR_AM5_PHYS.md).

**How to switch between branches?**
* `git checkout` allows you to switch between branches as long as they are present
in your local directory and any pending changes have been added via `git add` and `git commit`.
The checkout will change the files in your current directory to match the branch specified, so it is important to add and commit changes before using the command. If you have
a branch you would like to checkout, it will likely need to be fetched from the
remote repository it is in:
```bash
git fetch <remote repository name> <branch name>
git checkout <branch name>
```

**How to switch and merge branches during development?**
* Branches can also be used to store multiple copies of a branch without affecting existing developement on another branch.
If you were on a branch, addFeature, you can create a new branch with the existing history of your current branch:
```bash
## git add and commit any pending changes first
git checkout -b addFeature_testSomething
```
And then make some changes to test. You will then need to`git add` and `git commit` the updates to get
them onto the new branch. From there, you can test your changes to make sure they work
or if necessary `git push` your updates to store them in your forked repository.
Once you are done with the test branch, you can merge it into the original to add the
changes:
```bash
git checkout addFeature
git merge addFeature_testSomething
git push
```
When creating new branches, they will only be present locally on your system until the `git push <remote name> <branch name>` command is run.

#### Commits
**Best practices for commits**
* Commits should be small and have descriptive commit messages.  Refer to the
[Contributing guidelines](CONTRIBUTING.md) for more details on this.
* Please commit and push often. Small commits make differences easier to track. Pushing often allows for
easy sharing.

**How to view changes before adding them?**
* Before adding changes, `git diff` be used to show any current differences.

**How to get an older version from the commit history?**
* If you need an older version of the code, `git log` will show each commit added to the
repository in chronological order. To get a version of the code from the older commit
you would copy the hash (long random string after the word commit) and use it to checkout.
For example:
```bash
git checkout 1b5e9c2496596510a59bc1f6f7f1005721d36a4a
```


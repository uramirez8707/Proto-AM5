# Merging Code into this Repository

For FAQ and an overview on git please see the [FAQ document](FAQ.md).

### Set up fork and local repository

To contribute to this repository, a fork of this repository must be created under your account.
If unfamilar, this is a copy of the repository specific to your account to allow for changes without
affecting the original repository.

1. Create a fork of the AM5_phys repo (you will only do this once)
- Navigate to https://gitlab.gfdl.noaa.gov/fms/am5_phys in your web browser (log in using
GFDL credentials)
- Click on `fork` on the top right hand corner
- In the `Select a namespace` tab, select your name
- In the `Visibility level` tab, click on public (it will only be visible to those inside the firewall)
- Click on the blue `Fork project` button
- You now have your own copy of the AM5_phys repo located at https://gitlab.gfdl.noaa.gov/{gitlab_username}/am5_phys where {gitlab_username} is your gitlab username

2. Add your fork locally

If you followed the instructions in the quickstart you will already have a local copy of the
repository in your src directory. You can then enter the repository directory and add your fork as a remote repository to get access to it.
This remote will allow you to pull, push and modify branches in your forked repository created in step 1.
Here the added remote is given the name 'myfork', and by default the originally cloned
remote will already be present and given the name 'origin':

```
cd am5_phys
git remote add myfork https://gitlab.gfdl.noaa.gov/{gitlab_username}/am5_phys.git
```

If you haven't followed the quickstart, you can clone the original repository to get the most up-to-date
version, and then add the remote to access your fork:
```
git clone https://gitlab.gfdl.noaa.gov/fms/am5_phys.git
cd am5_phys
git add myfork https://gitlab.gfdl.noaa.gov/{gitlab_username}/am5_phys.git
```

After adding remotes, `git remote -v` can be used to check remote names and corresponding URL's.
For example, you should see something similar to:
```
myfork  https://gitlab.gfdl.noaa.gov/First.Last/am5_phys.git (fetch)
myfork	https://gitlab.gfdl.noaa.gov/First.Last/am5_phys.git (push)
origin  https://gitlab.gfdl.noaa.gov/m5/am5_phys.git (fetch)
origin	https://gitlab.gfdl.noaa.gov/m5/am5_phys.git (push)

```

3. Create a branch with a descriptive name

By default you will be on the `main` branch, but a new branch should be created and switched to
specifically for the changes you would like to add. To create and switch to a new local branch with
the name 'myUpdate':
```
git checkout -b myUpdate
```

For additional information on branching and remotes, please see the [best practices document](BEST_PRACTICES.md).

### Adding your changes

4. Make your updates
- You can look at the changes you made using
```
git diff
```

- You can look at the files that you changed using
```
git status
```

These will only display changes for files already present in the repository. Any newly created files
will not be displayed until added with the `git add` command below.

5. Add and commit your updates
```
git add <names or paths of all added and modified files>
git commit -m "Descriptive commit message saying exactly what was done"
```
You can view your commit in the branch's commit history with `git log`.

6. Make sure that your branch is up to date with the main branch
```
git fetch origin main
git merge origin/main
```
If there are merge conflicts, they will need to resolved. See the [Resolving merge conflicts](MERGE_FOR_XML.md#Resolving merge conflicts) guide.

7. Push updates to your own fork
```
git push myfork myUpdate
```
You will be asked login using the same GFDL username and pass used to log in to the gitlab site.

### Creating a Merge request

8. Navigate to https://gitlab.gfdl.noaa.gov/{gitlab_username}/am5_phys/-/merge_requests/new in your web browser, where {gitlab_username} is your gitlab username.

9. Select your newly created branch as the 'Source Branch', and leave the target branch as the
`main` branch for `fms/am5_phys` (or another branch if needed). Then select 'Compare branches and continue'

10. Fill out Merge Request Information
Create a descriptive title of what changes are being applied with the request, and then proceed to
fill out the template given in the description, including the checklist. If there are any changes to
data files, please update the documentation and explain below. Any users you would like to review
the changes can be assigned in the drop down box below the descriptions, or assigned if any actions
are required.
![image](docs_images/mr_4.png)
Before submitting, commits and changes can be viewed at the bottom of the page. When done, scroll
down and click the `Create Merge Request` button.

11. After submitting you should navigate to https://gitlab.gfdl.noaa.gov/fms/am5_phys/-/merge_requests.
Here, your merge request should be listed among all open merge requests for the repository. If not,
your merge request may be in your forked repository instead, and a new one will need to be created while ensuring the target repository is set to `fms/am5_phys`.

Once the merge request is submitted succesfully and the checklist actions listed in the description are complete, please assign the merge request to Uriel.Ramirez or Raymond.Menzel

# Resolving merge conflicts
At some point, you will receive messages that files have merge conflicts after pulling changes from the main branch. This happens when both branches/commits changed the same lines of code. These must be resolved before your code is merged to the main branch.

Git will show any conflicting areas between <<<<< and >>>>> with branches seperated by ======,
and include branch names for which branch the changes are from. You then need to pick which set of
code to keep or make any modifications to ensure the conflict is resolved.

- `git status` will show files that need to be resolved.
- Open the file with your favorite text editor
- Fix the conflicts and delete the <<<<< ====== >>>>>> lines
- Add and commit your changes
```
git add .
git commit -m "merge main and solved merge conflicts"
```

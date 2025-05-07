# AM5_Phys
This is the repository that contains AM5 physics!
## How to get started?
- ### [Quickstart for running AM5 models](https://gitlab.gfdl.noaa.gov/m5/am5xml#quickstart)
	-  Commands for running models, post-processing and testing with this code base with the AM5 XML
- ### [Quickstart for submitting AM5_physics **code** merge requests](MERGE_FOR_AM5_PHYS.md)
	- Please also see the repository [policies](POLICIES.md) and [contribution guide](CONTRIBUTING.md) before submitting a merge request
	- For updates to the am5 **xmls** please see the [instructions](https://gitlab.gfdl.noaa.gov/m5/am5xml/-/blob/main/MERGE_FOR_XML.md) for that repository.

## Documentation and Policies
- ### [Reproducibility Checking](https://gitlab.gfdl.noaa.gov/m5/am5xml/-/blob/main/REPRODUCIBILITY.md)
- ### [Contributing Guidelines](CONTRIBUTING.md)
- ### [Style Guide](STYLE.md)
- ### [FRE and XML tips](https://gitlab.gfdl.noaa.gov/m5/am5xml#fre-and-xml-tips-1)
- ### [Policies](POLICIES.md)
- ### Git and Gitlab documentation
	- ### [Git and Gitlab FAQ and Best Practices](FAQ.md)
	- ### [MSD General Guide to git](https://sites.google.com/noaa.gov/oar-gfdl-msd-docs/github-gitlab/git-guide-and-best-practices)
- ### [GFDL Fair Use Policy](#gfdl-fair-use-policy-for-code-and-data)
Please visit the above link for the entire fair use policy.  Below is the appendix related to Experimental Models

### GFDL Fair Use Policy for Experimental GFDL models
Experimental GFDL models refers to models (e.g., a coupled climate model), model components (e.g., parameterization schemes), and model configurations (e.g., the specific arrangement and parameter settings of model components) arising from model development efforts at GFDL, but whose formulation and configuration have not yet been documented in the peer-reviewed literature, as well as outputs from simulations with these models. Those wishing to make use of such data should contact the model developers to discuss potential model use, and request approval from the developers for planned use of the data. It is strongly desired that such use would be in the form of a collaboration with the developers. Any products derived from the model use (papers, presentations, etc.) should give appropriate credit to both the model developers and the collaborators.

## Using am5_phys in existing xmls
There are a few modifications to an XML in order switch to the am5_phys code. 
An XML will have a section in the compile experiment that looks like the following:
```xml
    <component name="atmos_phys" requires="fms" paths="atmos_phys">
      <description domainName="" communityName="" communityVersion="$(RELEASE)" communityGrid=""/>
      <source versionControl="git" root="http://gitlab.gfdl.noaa.gov/fms">
        <codeBase version="$(RELEASE)"> atmos_phys.git </codeBase>
          <csh><![CDATA[
            ( cd atmos_phys  && git checkout $(ATMOS_GIT_TAG) )
           ]]>
          </csh>
      </source>
      <compile>
        <cppDefs>$(F2003_FLAGS) -DCLUBB</cppDefs>
      </compile>
    </component>
```
To use AM5 physics you can do
```xml
    <component name="atmos_phys" requires="fms" paths="am5_phys">
      <description domainName="" communityName="" communityVersion="$(RELEASE)" communityGrid=""/>
        <source versionControl="git" root="http://gitlab.gfdl.noaa.gov/fms">
          <codeBase version="$(AM5_GIT_TAG)">am5_phys.git</codeBase>
        </source>
      <compile>
        <cppDefs>$(F2003_FLAGS) -DCLUBB </cppDefs>
      </compile>
    </component>
```
1. The **paths** has been changed to "am5_phys".
2. The **codeBase** was changed from **atmos_phys.git** to **am5_phys.git**
3. The `csh` block can be removed and checkout the tag you want by defining `AM5_GIT_TAG`

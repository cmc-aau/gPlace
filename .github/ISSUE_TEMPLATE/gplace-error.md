---
name: gPlace error
about: Help is required for a gPlace error.
title: ''
labels: error
assignees: ''

---

<!---
Resources for checking if the issue has already been addressed:
FAQ: https://cmc-aau.github.io/gPlace/faq.html
Issues: https://github.com/cmc-aau/gPlace/issues

Otherwise, please provide as much of the requested information below as possible.  
-->

## Environment
- [ ] Installed via pip (include the output of `pip list`)
- [ ] Using a conda environment (include the output of `conda list && conda list --revisions`)
- [ ] Using a Docker container (include the `IMAGE ID` of the container)

## Server information
- CPU: `grep -m 1 "^model name" /proc/cpuinfo && grep -c "^processor" /proc/cpuinfo`
- RAM: `grep "^MemTotal" /proc/meminfo`
- OS: `cat /etc/os-release`

## Debugging information
- [ ] `gplace.log` has been included (drag and drop the file to upload).
- [ ] MSA has been included (if possible).

## Additional comments
<!--- Include any additional comments which may help with debugging here. -->

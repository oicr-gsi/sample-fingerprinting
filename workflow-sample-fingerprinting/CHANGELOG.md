## 2.4.2 - 2019-08-19
- [GP-2060](https://jira.oicr.on.ca/browse/GP-2060?filter=-1)Fixed cluster.r script to address issue when wf does not produce plots
## 2.4.1 - 2018-11-01
- [GP-1867](https://jira.oicr.on.ca/browse/GP-1867?filter=-1)Fixed a regular expression in make_report.pl script preventing processing such studies as TGL08 or C4GD (numbers in donor's names)
## 2.4   - 2016-12-13
- Added two new jobs - optional alternative to jaccard matric building and proximity table generator script (specs for table in perldoc of make_table.pl script)
## 2.2.2 - 2016-01-19
- Introduced another solution to older tbi indices (tbi files are copied instead of linked to)
## 2.2.1 - 2015-09-30
- Added fix for removing leading zeros (was appearing with small number of samples that was possible to fit on one heatmap)
## 2.2   - 2015-09-09
- Added additional code for attaching annotations to result files showing which set of hotspots to use
## 2.1.2 - 2015-08-14
- Fixed an issue with executing touch, adjusted memory for matrix-producing step, regex for extracting donor name
## 2.1.1 - 2015-07-31
- Had to fix a regex in new swap detection code that could miss some branches in dendrogram b/c of scientific notation for the distance value
## 2.1   - 2015-07-24
- New Swap-detection algorithm based on dendrogram analysis
## 2.0.2 - 2015-05-25
- Fixed paths for test data
## 2.0.1 - 2015-03-19
- Upgrade to SeqWare 1.1.0, common-utilities 1.6 and workflow-utilities 1.6.

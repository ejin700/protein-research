This program datamines the UniProt website for sequences that have 
experimentally proven n-linked glycosylation sites. It can also perform some statistical analysis on the output data. To use this program:

1). Run info_tester.py
	-requires BioPython, tqdm, Bioservices packages
	-SPECIFY PATH!
	-can alter inital query or evidence terms

2). Run alignment_tester.py
	-make sure clustalo and extentions are in working directory
	-SPECIFY PATH! (same as above)

3). Run ansifilter from command line on text files in 'pretty' folder
	-first, copy all files to another folder
	-make sure ansifiletr and extentions in working dircetory
	-command is 'ansifilter *.txt -H'
	-delete all text files (keep all html files)
	-to rename html files, command is 'Dir | Rename-Item -NewName {$_.name -replace ".txt",""}

For the statistical analysis (4), I initially wrote each function in a different file, then copy and pasted them into analysis.py and did my best to convert it into a module file. Also, now that all the functions are together in analysis.py, we can definitely write new functions that can be used as intermediates in the other functions. The analysis_tester.py is not really a runable files; it just has separate function calls from analysis.py and a very rough pipeline structure. 
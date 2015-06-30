##AuAu14.6GeV Omg Embedding Notes
The full statistics of Omega and AntiOmega and Phi mesons for AuAu 14.5GeV embedding data is ready.       
-20150624;

###How to submit jobs on pdsf
1. Need to modify the makelist.py to generate filelist of full statistics;
1. You need to update your `StRefMultCorr/` directory under your `StRoot`(how to check out the most updated one?);
1. Change the GeantId for omega and lambda;(why lambda has two different geant ids)
1. Replace my StMcV0Maker and StMcXiMaker with Feng's
1. Change the `StRoot/StV0Maker/StMcV0Maker.cxx`'s trigger ids;
1. Change the starver by `$starver SL14i`;
1. Change the file list name in `qsub_script_omg.sh`;
1. Change `run.csh`, `run2.csh`, `submit_omg.sh`, `qsub_script_omg.sh` working directory;
1. Revise `qsub_script_omg.sh` to

```bash
#!/bin/bash
cd /global/homes/l/lwen1990/pwg/embedding/15GeV_embedding
qsub -t 1-70:1 submit_omg.sh #the number 70 can be calculated by ceil(No. of total files / 10)
```

  i. to submit jobs, use `./qsub_script_omg.sh`
  
  ii. you can use `qstat -u lwen1990` to check your jobs status
  
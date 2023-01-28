
# amplicon_sorter

THIS IS A DRAFT VERSION AND NOT THOROUGHLY TESTED.  USE AT OWN RISK.

### Extra Added Options:

`-ar', '--allreads`: Use all available reads between the length limits to process.  This argument is still limited with "-maxreads" to have a hard limit.

`-ldc', '--length_diff_consensus`: Length difference between consensusses allowed to COMBINE groups based on the consensus sequence (value between 0 and 100). Default=8.0%.  It is possible that you will have more than 1 group with the same sequence, but different in length (caused by incomplete PCR).  By default Amplicon_sorter will keep these groups separate.  If you want those groups merged for some reason, you can use this option.  The consensus will probably be shorter (the average of the several lengths).  Still in testing phase.




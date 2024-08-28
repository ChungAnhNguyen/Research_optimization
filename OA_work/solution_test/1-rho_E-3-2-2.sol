<?xml version = "1.0" encoding="UTF-8" standalone="yes"?>
<CPLEXSolution version="1.2">
 <header
   problemName="1_rho__E_3_2_2"
   solutionName="incumbent"
   solutionIndex="-1"
   objectiveValue="4"
   solutionTypeValue="3"
   solutionTypeString="primal"
   solutionStatusValue="101"
   solutionStatusString="integer optimal solution"
   solutionMethodString="mip"
   primalFeasible="1"
   dualFeasible="1"
   MIPNodes="0"
   MIPIterations="1"
   writeLevel="1"/>
 <quality
   epInt="1.0000000000000001e-05"
   epRHS="9.9999999999999995e-07"
   maxIntInfeas="0"
   maxPrimalInfeas="0"
   maxX="4"
   maxSlack="0"/>
 <linearConstraints>
  <constraint name="R_def" index="0" slack="0"/>
  <constraint name="val_cible" index="1" slack="0"/>
  <constraint name="BkI_6_0" index="2" slack="0"/>
  <constraint name="BkI_6_1" index="3" slack="0"/>
  <constraint name="BkI_5_0" index="4" slack="0"/>
  <constraint name="BkI_5_1" index="5" slack="0"/>
  <constraint name="BkI_3_0" index="6" slack="0"/>
  <constraint name="BkI_3_1" index="7" slack="0"/>
  <constraint name="QstargeqQ1" index="8" slack="0"/>
  <constraint name="QstargeqQ2" index="9" slack="0"/>
  <constraint name="QstargeqQ3" index="10" slack="0"/>
 </linearConstraints>
 <variables>
  <variable name="R" index="0" value="4"/>
  <variable name="Q000" index="1" value="1"/>
  <variable name="Q001" index="2" value="1"/>
  <variable name="Q010" index="3" value="1"/>
  <variable name="Q011" index="4" value="1"/>
 </variables>
 <objectiveValues>
  <objective index="0" name="" value="4"/>
 </objectiveValues>
</CPLEXSolution>

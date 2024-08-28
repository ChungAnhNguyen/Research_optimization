<?xml version = "1.0" encoding="UTF-8" standalone="yes"?>
<CPLEXSolution version="1.2">
 <header
   problemName="0_gamma__E_3_2_2"
   objectiveValue="0"
   solutionTypeValue="1"
   solutionTypeString="basic"
   solutionStatusValue="3"
   solutionStatusString="infeasible"
   solutionMethodString="dual"
   primalFeasible="0"
   dualFeasible="1"
   simplexIterations="2"
   writeLevel="1"/>
 <quality
   epRHS="9.9999999999999995e-07"
   epOpt="9.9999999999999995e-07"
   maxPrimalInfeas="2"
   maxDualInfeas="0"
   maxPrimalResidual="0"
   maxDualResidual="5.5511151231257827e-17"
   maxX="1"
   maxPi="0.33333333333333331"
   maxSlack="2"
   maxRedCost="0"
   kappa="80.25"/>
 <linearConstraints>
  <constraint name="R_def" index="0" status="BS" slack="0" dual="-0"/>
  <constraint name="val_cible_G" index="1" status="BS" slack="2" dual="0"/>
  <constraint name="BkI_3_0" index="2" status="BS" slack="0" dual="-0"/>
  <constraint name="BkI_3_1" index="3" status="LL" slack="0" dual="-0.33333333333333331"/>
  <constraint name="BkI_3_2" index="4" status="LL" slack="0" dual="0.33333333333333331"/>
  <constraint name="BkI_5_0" index="5" status="BS" slack="0" dual="-0"/>
  <constraint name="BkI_5_1" index="6" status="LL" slack="0" dual="0.33333333333333331"/>
  <constraint name="BkI_5_2" index="7" status="LL" slack="0" dual="-0.33333333333333331"/>
  <constraint name="BkI_6_0" index="8" status="BS" slack="0" dual="-0"/>
  <constraint name="BkI_6_1" index="9" status="LL" slack="0" dual="-0.33333333333333331"/>
  <constraint name="BkI_6_2" index="10" status="LL" slack="0" dual="0.33333333333333331"/>
 </linearConstraints>
 <variables>
  <variable name="R" index="0" status="LL" value="1" reducedCost="-0"/>
  <variable name="Q012" index="1" status="BS" value="0" reducedCost="0"/>
  <variable name="Q000" index="2" status="UL" value="1" reducedCost="0"/>
  <variable name="P000" index="3" status="UL" value="1" reducedCost="0"/>
  <variable name="Q001" index="4" status="LL" value="0" reducedCost="-0"/>
  <variable name="P001" index="5" status="BS" value="0" reducedCost="0"/>
  <variable name="Q002" index="6" status="BS" value="0" reducedCost="0"/>
  <variable name="P002" index="7" status="LL" value="0" reducedCost="-0"/>
  <variable name="Q010" index="8" status="BS" value="0" reducedCost="0"/>
  <variable name="P010" index="9" status="LL" value="0" reducedCost="-0"/>
  <variable name="Q011" index="10" status="BS" value="0" reducedCost="0"/>
  <variable name="P011" index="11" status="LL" value="0" reducedCost="-0"/>
  <variable name="Q020" index="12" status="LL" value="0" reducedCost="-0"/>
  <variable name="P020" index="13" status="LL" value="0" reducedCost="-0"/>
  <variable name="Q022" index="14" status="BS" value="0" reducedCost="0"/>
  <variable name="P022" index="15" status="LL" value="0" reducedCost="-0"/>
 </variables>
 <objectiveValues>
  <objective index="0" name="" value="0"/>
 </objectiveValues>
</CPLEXSolution>

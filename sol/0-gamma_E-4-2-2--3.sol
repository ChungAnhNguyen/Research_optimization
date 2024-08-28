<?xml version = "1.0" encoding="UTF-8" standalone="yes"?>
<CPLEXSolution version="1.2">
 <header
   problemName="0_gamma__E_4_2_2"
   objectiveValue="0.25000000000000006"
   solutionTypeValue="1"
   solutionTypeString="basic"
   solutionStatusValue="1"
   solutionStatusString="optimal"
   solutionMethodString="dual"
   primalFeasible="1"
   dualFeasible="1"
   simplexIterations="4"
   writeLevel="1"/>
 <quality
   epRHS="9.9999999999999995e-07"
   epOpt="9.9999999999999995e-07"
   maxPrimalInfeas="1.1102230246251565e-16"
   maxDualInfeas="6.6613381477509392e-16"
   maxPrimalResidual="1.1102230246251565e-16"
   maxDualResidual="7.2164496600635175e-16"
   maxX="1"
   maxPi="0.75000000000000022"
   maxSlack="2.4999999999999991"
   maxRedCost="0.75000000000000022"
   kappa="41.437773722627725"/>
 <linearConstraints>
  <constraint name="R_def" index="0" status="LL" slack="0" dual="0.25000000000000006"/>
  <constraint name="val_cible_G" index="1" status="BS" slack="-2.4999999999999991" dual="0"/>
  <constraint name="BkI_3_0" index="2" status="BS" slack="8.3266726846886741e-17" dual="-0"/>
  <constraint name="BkI_3_1" index="3" status="LL" slack="0" dual="-0.25000000000000006"/>
  <constraint name="BkI_3_2" index="4" status="LL" slack="0" dual="-1.6653345369377348e-16"/>
  <constraint name="BkI_3_3" index="5" status="LL" slack="0" dual="2.2204460492503131e-16"/>
  <constraint name="BkI_5_0" index="6" status="LL" slack="0" dual="0.062499999999999972"/>
  <constraint name="BkI_5_1" index="7" status="BS" slack="-2.7755575615628914e-17" dual="-0"/>
  <constraint name="BkI_5_2" index="8" status="LL" slack="0" dual="-0.062499999999999944"/>
  <constraint name="BkI_5_3" index="9" status="LL" slack="0" dual="0.12499999999999983"/>
  <constraint name="BkI_6_0" index="10" status="LL" slack="0" dual="0.062500000000000083"/>
  <constraint name="BkI_6_1" index="11" status="BS" slack="-2.7755575615628914e-17" dual="-0"/>
  <constraint name="BkI_6_2" index="12" status="LL" slack="0" dual="0.18750000000000028"/>
  <constraint name="BkI_6_3" index="13" status="LL" slack="0" dual="0.12500000000000017"/>
  <constraint name="BkI_9_0" index="14" status="LL" slack="0" dual="0.31249999999999994"/>
  <constraint name="BkI_9_1" index="15" status="LL" slack="0" dual="0.37500000000000011"/>
  <constraint name="BkI_9_2" index="16" status="LL" slack="0" dual="0.18750000000000061"/>
  <constraint name="BkI_9_3" index="17" status="BS" slack="-2.7755575615628914e-17" dual="-0"/>
  <constraint name="BkI_10_0" index="18" status="LL" slack="0" dual="-0.5625"/>
  <constraint name="BkI_10_1" index="19" status="LL" slack="0" dual="-0.75000000000000022"/>
  <constraint name="BkI_10_2" index="20" status="LL" slack="0" dual="-0.6875"/>
  <constraint name="BkI_10_3" index="21" status="LL" slack="0" dual="-0.62500000000000011"/>
  <constraint name="BkI_12_0" index="22" status="LL" slack="0" dual="0.125"/>
  <constraint name="BkI_12_1" index="23" status="BS" slack="-1.1102230246251565e-16" dual="-0"/>
  <constraint name="BkI_12_2" index="24" status="LL" slack="0" dual="0.12500000000000003"/>
  <constraint name="BkI_12_3" index="25" status="LL" slack="0" dual="0.24999999999999972"/>
 </linearConstraints>
 <variables>
  <variable name="R" index="0" status="UL" value="1" reducedCost="0.25000000000000006"/>
  <variable name="Q0123" index="1" status="BS" value="0.24999999999999994" reducedCost="0"/>
  <variable name="Q0000" index="2" status="BS" value="0.25000000000000006" reducedCost="0"/>
  <variable name="P0000" index="3" status="LL" value="0" reducedCost="-0.25"/>
  <variable name="Q0001" index="4" status="LL" value="0" reducedCost="-0.24999999999999989"/>
  <variable name="P0001" index="5" status="BS" value="0" reducedCost="0"/>
  <variable name="Q0002" index="6" status="LL" value="0" reducedCost="-0.24999999999999944"/>
  <variable name="P0002" index="7" status="LL" value="0" reducedCost="-6.106226635438361e-16"/>
  <variable name="Q0003" index="8" status="LL" value="0" reducedCost="-0.24999999999999989"/>
  <variable name="P0003" index="9" status="LL" value="0" reducedCost="-2.2204460492503131e-16"/>
  <variable name="Q0010" index="10" status="LL" value="0" reducedCost="-0.25000000000000011"/>
  <variable name="P0010" index="11" status="LL" value="0" reducedCost="-0"/>
  <variable name="Q0011" index="12" status="LL" value="0" reducedCost="-0.25000000000000011"/>
  <variable name="P0011" index="13" status="BS" value="0.24999999999999997" reducedCost="0"/>
  <variable name="Q0012" index="14" status="LL" value="0" reducedCost="-0.74999999999999956"/>
  <variable name="Q0013" index="15" status="LL" value="0" reducedCost="-0.50000000000000056"/>
  <variable name="Q0020" index="16" status="LL" value="0" reducedCost="-0"/>
  <variable name="P0020" index="17" status="LL" value="0" reducedCost="-0.25000000000000022"/>
  <variable name="Q0021" index="18" status="LL" value="0" reducedCost="6.6613381477509392e-16"/>
  <variable name="Q0022" index="19" status="LL" value="0" reducedCost="-0.24999999999999906"/>
  <variable name="P0022" index="20" status="LL" value="0" reducedCost="-9.9920072216264089e-16"/>
  <variable name="Q0023" index="21" status="LL" value="0" reducedCost="-0.5"/>
  <variable name="Q0030" index="22" status="LL" value="0" reducedCost="-0.25000000000000011"/>
  <variable name="P0030" index="23" status="BS" value="0" reducedCost="0"/>
  <variable name="Q0031" index="24" status="BS" value="0.24999999999999994" reducedCost="0"/>
  <variable name="Q0032" index="25" status="LL" value="0" reducedCost="-0.24999999999999906"/>
  <variable name="Q0033" index="26" status="LL" value="0" reducedCost="-0.25000000000000011"/>
  <variable name="P0033" index="27" status="BS" value="0.24999999999999997" reducedCost="0"/>
  <variable name="Q0100" index="28" status="LL" value="0" reducedCost="-0"/>
  <variable name="P0100" index="29" status="LL" value="0" reducedCost="-0.25"/>
  <variable name="Q0101" index="30" status="LL" value="0" reducedCost="-0.24999999999999978"/>
  <variable name="P0101" index="31" status="BS" value="0.25000000000000006" reducedCost="0"/>
  <variable name="Q0102" index="32" status="LL" value="0" reducedCost="-0.24999999999999939"/>
  <variable name="Q0103" index="33" status="LL" value="0" reducedCost="-0.49999999999999956"/>
  <variable name="Q0110" index="34" status="LL" value="0" reducedCost="-0.25000000000000011"/>
  <variable name="P0110" index="35" status="BS" value="0" reducedCost="0"/>
  <variable name="Q0111" index="36" status="LL" value="0" reducedCost="-0.25"/>
  <variable name="P0111" index="37" status="LL" value="0" reducedCost="-0"/>
  <variable name="Q0112" index="38" status="LL" value="0" reducedCost="-0.74999999999999956"/>
  <variable name="Q0113" index="39" status="LL" value="0" reducedCost="-0.75000000000000022"/>
  <variable name="Q0120" index="40" status="LL" value="0" reducedCost="-0.25000000000000039"/>
  <variable name="Q0121" index="41" status="LL" value="0" reducedCost="-0.24999999999999967"/>
  <variable name="Q0122" index="42" status="LL" value="0" reducedCost="-0.49999999999999939"/>
  <variable name="Q0130" index="43" status="LL" value="0" reducedCost="-0.25000000000000011"/>
  <variable name="Q0131" index="44" status="LL" value="0" reducedCost="1.1102230246251565e-16"/>
  <variable name="Q0133" index="45" status="LL" value="0" reducedCost="-0.49999999999999978"/>
  <variable name="Q0200" index="46" status="BS" value="0" reducedCost="0"/>
  <variable name="P0200" index="47" status="LL" value="0" reducedCost="-0.25000000000000033"/>
  <variable name="Q0201" index="48" status="LL" value="0" reducedCost="1.9428902930940239e-16"/>
  <variable name="Q0202" index="49" status="BS" value="0" reducedCost="0"/>
  <variable name="P0202" index="50" status="LL" value="0" reducedCost="-0.25000000000000089"/>
  <variable name="Q0203" index="51" status="LL" value="0" reducedCost="-0.24999999999999947"/>
  <variable name="Q0210" index="52" status="LL" value="0" reducedCost="-0.24999999999999981"/>
  <variable name="Q0211" index="53" status="BS" value="0.24999999999999997" reducedCost="0"/>
  <variable name="Q0212" index="54" status="LL" value="0" reducedCost="-0.49999999999999956"/>
  <variable name="Q0220" index="55" status="LL" value="0" reducedCost="-0.25000000000000011"/>
  <variable name="P0220" index="56" status="BS" value="0.24999999999999994" reducedCost="0"/>
  <variable name="Q0221" index="57" status="BS" value="0" reducedCost="0"/>
  <variable name="Q0222" index="58" status="LL" value="0" reducedCost="-0.24999999999999931"/>
  <variable name="P0222" index="59" status="LL" value="0" reducedCost="-6.9388939039072284e-16"/>
  <variable name="Q0223" index="60" status="LL" value="0" reducedCost="-0.75"/>
  <variable name="Q0230" index="61" status="LL" value="0" reducedCost="-0.50000000000000022"/>
  <variable name="Q0232" index="62" status="LL" value="0" reducedCost="-0.24999999999999931"/>
  <variable name="Q0233" index="63" status="LL" value="0" reducedCost="-0.50000000000000011"/>
  <variable name="Q0300" index="64" status="LL" value="0" reducedCost="-3.3306690738754696e-16"/>
  <variable name="P0300" index="65" status="LL" value="0" reducedCost="-0.24999999999999972"/>
  <variable name="Q0301" index="66" status="LL" value="0" reducedCost="-0.25000000000000017"/>
  <variable name="Q0302" index="67" status="BS" value="0" reducedCost="0"/>
  <variable name="Q0303" index="68" status="LL" value="0" reducedCost="-0.25000000000000006"/>
  <variable name="P0303" index="69" status="BS" value="0" reducedCost="0"/>
  <variable name="Q0310" index="70" status="BS" value="0" reducedCost="0"/>
  <variable name="Q0311" index="71" status="BS" value="0" reducedCost="0"/>
  <variable name="Q0313" index="72" status="LL" value="0" reducedCost="-0.25000000000000039"/>
  <variable name="Q0320" index="73" status="LL" value="0" reducedCost="-3.3306690738754696e-16"/>
  <variable name="Q0322" index="74" status="LL" value="0" reducedCost="3.3306690738754696e-16"/>
  <variable name="Q0323" index="75" status="LL" value="0" reducedCost="-0.50000000000000022"/>
  <variable name="Q0330" index="76" status="LL" value="0" reducedCost="-0.25000000000000039"/>
  <variable name="P0330" index="77" status="BS" value="0" reducedCost="0"/>
  <variable name="Q0331" index="78" status="LL" value="0" reducedCost="-2.7755575615628914e-16"/>
  <variable name="Q0332" index="79" status="BS" value="0" reducedCost="0"/>
  <variable name="Q0333" index="80" status="LL" value="0" reducedCost="-0.25000000000000028"/>
  <variable name="P0333" index="81" status="LL" value="0" reducedCost="2.7755575615628914e-16"/>
 </variables>
 <objectiveValues>
  <objective index="0" name="" value="0.24999999999999994"/>
 </objectiveValues>
</CPLEXSolution>

"""
Chapter 6, Hands-on 5 - Bayesian dating with BEAST 2: BETS and Bayes factors for clock models on a stratified subset of the Ebola virus genomes of Chapter 5 (by default 60 of 195, evenly spaced by collection date).
Writes BEAST 2.7 XML files (the standard BEAUti structure) for these models, all with HKY + 4 gamma categories (empirical base frequencies) and an exponential-growth coalescent tree prior:
  het_strict      real tip dates, strict clock                      (the model of interest)
  het_relaxed     real tip dates, uncorrelated lognormal relaxed clock
  iso_strict      no dates: all tips given the same date (isochronous)
  dr1, dr2, dr3   dates randomly permuted among the tips (date randomisation, strict clock)
Each model is run twice: as an ordinary MCMC (xml/mcmc_<model>.xml) for the posterior, and wrapped in a path sampler (xml/ps_<model>.xml, MODEL_SELECTION package) for the marginal likelihood by path sampling and stepping-stone sampling.
Priors are identical across models: clock rate LogNormal(M = -6.9, S = 1.25) per site per year; kappa LogNormal(1, 1.25); gamma shape Exponential(1); population size at the present LogNormal(M = 0, S = 3) years; growth rate Laplace(0, 30.7) per year.
USAGE   python ex6_beast.py make [--mcmc-length N] [--ps-steps K] [--ps-length N] [--ps-pre N] [--taxa N]      writes xml/*.xml
Needs: BEAST 2.7.7 on the path (beast) with the MODEL_SELECTION package 1.6.2, numpy, pandas.
"""
import argparse, datetime, os, sys
import numpy as np, pandas as pd
ap = argparse.ArgumentParser(); ap.add_argument("cmd"); ap.add_argument("--mcmc-length", type=int, default=20000000); ap.add_argument("--ps-steps", type=int, default=24); ap.add_argument("--ps-length", type=int, default=400000); ap.add_argument("--ps-pre", type=int, default=200000); ap.add_argument("--taxa", type=int, default=60); A = ap.parse_args()
def readfa(p):
    n, s = [], []
    for l in open(p):
        l = l.strip()
        if l.startswith(">"): n.append(l[1:].split()[0]); s.append([])
        elif l: s[-1].append(l.upper())
    return n, ["".join(x) for x in s]
names, seqs = readfa("ebola_aln.fasta"); dts = pd.read_csv("ebola_dates.csv").set_index("acc").date.to_dict()
order = sorted(range(len(names)), key=lambda i: (dts[names[i]], names[i])); pick = sorted(order[int(round(k))] for k in np.linspace(0, len(order) - 1, A.taxa)); names, seqs = [names[i] for i in pick], [seqs[i] for i in pick]      # every k-th genome by collection date
def dec(d): d = pd.Timestamp(d); return d.year + (d.dayofyear - 0.5) / (366 if d.is_leap_year else 365)
real = {a: dec(dts[a]) for a in names}; first = min(real.values()); last = max(real.values())
def dates_for(kind):
    if kind.startswith("het"): return real
    if kind == "iso": return {a: last for a in names}
    k = int(kind[2:]); rng = np.random.default_rng(100 + k); v = list(real.values()); p = rng.permutation(len(v)); return {a: v[p[i]] for i, a in enumerate(names)}
NS = "beast.pkgmgmt:beast.base.core:beast.base.inference:beast.base.evolution.alignment:beast.base.evolution.tree.coalescent:beast.base.util:beast.base.math:beast.evolution.nuc:beast.base.evolution.operator:beast.base.inference.operator:beast.base.evolution.sitemodel:beast.base.evolution.substitutionmodel:beast.base.evolution.likelihood"
def make(kind, mode, out):
    relaxed = kind == "het_relaxed"; dt = dates_for(kind); n = len(names); L = A.mcmc_length
    seq = "\n".join(f'        <sequence id="seq_{a}" spec="Sequence" taxon="{a}" totalcount="4" value="{s}"/>' for a, s in zip(names, seqs))
    trait = ",".join(f"{a}={dt[a]:.6f}" for a in names)
    if relaxed:
        clock = ('<branchRateModel id="RelaxedClock.c:ebola" spec="beast.base.evolution.branchratemodel.UCRelaxedClockModel" rateCategories="@rateCategories.c:ebola" tree="@Tree.t:ebola" clock.rate="@ucldMean.c:ebola">\n'
                 '                    <LogNormal id="LogNormalDistributionModel.c:ebola" S="@ucldStdev.c:ebola" meanInRealSpace="true" name="distr"><parameter id="RealParameter.M" spec="parameter.RealParameter" estimate="false" lower="0.0" name="M" upper="1.0" value="1.0"/></LogNormal>\n                </branchRateModel>')
        state_clock = f'        <stateNode id="rateCategories.c:ebola" spec="parameter.IntegerParameter" dimension="{2 * n - 2}" value="1"/>\n        <parameter id="ucldMean.c:ebola" spec="parameter.RealParameter" name="stateNode" value="0.001"/>\n        <parameter id="ucldStdev.c:ebola" spec="parameter.RealParameter" lower="0.0" name="stateNode" upper="5.0" value="0.3"/>'
        prior_clock = ('<prior id="ClockPrior.c:ebola" name="distribution" x="@ucldMean.c:ebola"><LogNormal id="LogNormalClock" name="distr" M="-6.9" S="1.25" meanInRealSpace="false"/></prior>\n'
                       '            <prior id="UcldStdevPrior.c:ebola" name="distribution" x="@ucldStdev.c:ebola"><Gamma id="GammaStdev" name="distr" alpha="0.5396" beta="0.3819"/></prior>')
        ops_clock = ('<operator id="ucldMeanScaler.c:ebola" spec="ScaleOperator" parameter="@ucldMean.c:ebola" scaleFactor="0.5" weight="2.0"/>\n    <operator id="ucldStdevScaler.c:ebola" spec="ScaleOperator" parameter="@ucldStdev.c:ebola" scaleFactor="0.5" weight="3.0"/>\n'
                     '    <operator id="CategoriesRandomWalk.c:ebola" spec="IntRandomWalkOperator" parameter="@rateCategories.c:ebola" weight="10.0" windowSize="1"/>\n    <operator id="CategoriesSwapOperator.c:ebola" spec="SwapOperator" intparameter="@rateCategories.c:ebola" weight="10.0"/>\n    <operator id="CategoriesUniform.c:ebola" spec="UniformOperator" parameter="@rateCategories.c:ebola" weight="10.0"/>\n'
                     '    <operator id="relaxedUpDownOperator.c:ebola" spec="UpDownOperator" scaleFactor="0.75" weight="3.0"><up idref="ucldMean.c:ebola"/><down idref="Tree.t:ebola"/></operator>')
        log_clock = '<log idref="ucldMean.c:ebola"/><log idref="ucldStdev.c:ebola"/><log id="rate.c:ebola" spec="beast.base.evolution.RateStatistic" branchratemodel="@RelaxedClock.c:ebola" tree="@Tree.t:ebola"/>'; rate_name = "ucldMean.c:ebola"
    else:
        clock = '<branchRateModel id="StrictClock.c:ebola" spec="beast.base.evolution.branchratemodel.StrictClockModel"><parameter id="clockRate.c:ebola" spec="parameter.RealParameter" estimate="true" lower="0.0" name="clock.rate" value="0.001"/></branchRateModel>'
        state_clock = '        <parameter id="clockRate.c:ebola" spec="parameter.RealParameter" lower="0.0" name="stateNode" value="0.001"/>'
        clock = '<branchRateModel id="StrictClock.c:ebola" spec="beast.base.evolution.branchratemodel.StrictClockModel" clock.rate="@clockRate.c:ebola"/>'
        prior_clock = '<prior id="ClockPrior.c:ebola" name="distribution" x="@clockRate.c:ebola"><LogNormal id="LogNormalClock" name="distr" M="-6.9" S="1.25" meanInRealSpace="false"/></prior>'
        ops_clock = '<operator id="StrictClockRateScaler.c:ebola" spec="ScaleOperator" parameter="@clockRate.c:ebola" scaleFactor="0.75" weight="3.0"/>\n    <operator id="strictClockUpDownOperator.c:ebola" spec="UpDownOperator" scaleFactor="0.75" weight="3.0"><up idref="clockRate.c:ebola"/><down idref="Tree.t:ebola"/></operator>'
        log_clock = '<log idref="clockRate.c:ebola"/>'
    mcmc_attr = f'chainLength="{L}"' if mode == "mcmc" else 'chainLength="1000000"'
    mcmc = f'''  <mcmc id="mcmc" spec="MCMC" {mcmc_attr}>
    <state id="state" spec="State" storeEvery="5000">
        <tree id="Tree.t:ebola" spec="beast.base.evolution.tree.Tree" name="stateNode">
            <trait id="dateTrait.t:ebola" spec="beast.base.evolution.tree.TraitSet" traitname="date" value="{trait}">
                <taxa id="TaxonSet.ebola" spec="TaxonSet"><alignment idref="ebola"/></taxa>
            </trait>
            <taxonset idref="TaxonSet.ebola"/>
        </tree>
        <parameter id="popSize.t:ebola" spec="parameter.RealParameter" name="stateNode" value="0.5"/>
        <parameter id="growthRate.t:ebola" spec="parameter.RealParameter" name="stateNode" value="10.0"/>
        <parameter id="kappa.s:ebola" spec="parameter.RealParameter" lower="0.0" name="stateNode" value="2.0"/>
        <parameter id="gammaShape.s:ebola" spec="parameter.RealParameter" name="stateNode" value="1.0"/>
{state_clock}
    </state>
    <init id="RandomTree.t:ebola" spec="beast.base.evolution.tree.coalescent.RandomTree" estimate="false" initial="@Tree.t:ebola" taxa="@ebola">
        <populationModel id="ConstantPopulation0.t:ebola" spec="ConstantPopulation"><parameter id="randomPopSize.t:ebola" spec="parameter.RealParameter" name="popSize" value="1.0"/></populationModel>
    </init>
    <distribution id="posterior" spec="beast.base.inference.CompoundDistribution">
        <distribution id="prior" spec="beast.base.inference.CompoundDistribution">
            <distribution id="CoalescentExponential.t:ebola" spec="Coalescent">
                <populationModel id="ExponentialGrowth.t:ebola" spec="ExponentialGrowth" growthRate="@growthRate.t:ebola" popSize="@popSize.t:ebola"/>
                <treeIntervals id="TreeIntervals.t:ebola" spec="beast.base.evolution.tree.TreeIntervals" tree="@Tree.t:ebola"/>
            </distribution>
            <prior id="PopSizePrior.t:ebola" name="distribution" x="@popSize.t:ebola"><LogNormal id="LogNormalPop" name="distr" M="0.0" S="3.0" meanInRealSpace="false"/></prior>
            <prior id="GrowthRatePrior.t:ebola" name="distribution" x="@growthRate.t:ebola"><LaplaceDistribution id="LaplaceGrowth" name="distr" mu="0.0" scale="30.701"/></prior>
            <prior id="KappaPrior.s:ebola" name="distribution" x="@kappa.s:ebola"><LogNormal id="LogNormalKappa" name="distr" M="1.0" S="1.25" meanInRealSpace="false"/></prior>
            <prior id="GammaShapePrior.s:ebola" name="distribution" x="@gammaShape.s:ebola"><Exponential id="ExpShape" name="distr" mean="1.0"/></prior>
            {prior_clock}
        </distribution>
        <distribution id="likelihood" spec="beast.base.inference.CompoundDistribution" useThreads="false">
            <distribution id="treeLikelihood.ebola" spec="TreeLikelihood" data="@ebola" tree="@Tree.t:ebola">
                <siteModel id="SiteModel.s:ebola" spec="SiteModel" gammaCategoryCount="4" shape="@gammaShape.s:ebola">
                    <parameter id="mutationRate.s:ebola" spec="parameter.RealParameter" estimate="false" name="mutationRate" value="1.0"/>
                    <parameter id="proportionInvariant.s:ebola" spec="parameter.RealParameter" estimate="false" lower="0.0" name="proportionInvariant" upper="1.0" value="0.0"/>
                    <substModel id="hky.s:ebola" spec="HKY" kappa="@kappa.s:ebola"><frequencies id="empiricalFreqs.s:ebola" spec="Frequencies" data="@ebola"/></substModel>
                </siteModel>
                {clock}
            </distribution>
        </distribution>
    </distribution>
    <operator id="CoalescentExponentialTreeScaler.t:ebola" spec="ScaleOperator" scaleFactor="0.5" tree="@Tree.t:ebola" weight="3.0"/>
    <operator id="CoalescentExponentialTreeRootScaler.t:ebola" spec="ScaleOperator" rootOnly="true" scaleFactor="0.5" tree="@Tree.t:ebola" weight="3.0"/>
    <operator id="CoalescentExponentialUniformOperator.t:ebola" spec="Uniform" tree="@Tree.t:ebola" weight="30.0"/>
    <operator id="CoalescentExponentialSubtreeSlide.t:ebola" spec="SubtreeSlide" tree="@Tree.t:ebola" weight="15.0"/>
    <operator id="CoalescentExponentialNarrow.t:ebola" spec="Exchange" tree="@Tree.t:ebola" weight="15.0"/>
    <operator id="CoalescentExponentialWide.t:ebola" spec="Exchange" isNarrow="false" tree="@Tree.t:ebola" weight="3.0"/>
    <operator id="CoalescentExponentialWilsonBalding.t:ebola" spec="WilsonBalding" tree="@Tree.t:ebola" weight="3.0"/>
    <operator id="PopSizeScaler.t:ebola" spec="ScaleOperator" parameter="@popSize.t:ebola" scaleFactor="0.75" weight="3.0"/>
    <operator id="GrowthRateRandomWalk.t:ebola" spec="RealRandomWalkOperator" parameter="@growthRate.t:ebola" weight="3.0" windowSize="1.0"/>
    <operator id="KappaScaler.s:ebola" spec="ScaleOperator" parameter="@kappa.s:ebola" scaleFactor="0.5" weight="0.1"/>
    <operator id="gammaShapeScaler.s:ebola" spec="ScaleOperator" parameter="@gammaShape.s:ebola" scaleFactor="0.5" weight="0.1"/>
    {ops_clock}
    <logger id="tracelog" spec="Logger" fileName="{os.path.basename(out)[:-4]}.log" logEvery="{max(L // 2000, 1000) if mode == 'mcmc' else 1000}" model="@posterior" sanitiseHeaders="true" sort="smart">
        <log idref="posterior"/><log idref="likelihood"/><log idref="prior"/><log idref="treeLikelihood.ebola"/>
        <log id="TreeHeight.t:ebola" spec="beast.base.evolution.tree.TreeHeightLogger" tree="@Tree.t:ebola"/>
        <log idref="popSize.t:ebola"/><log idref="growthRate.t:ebola"/><log idref="kappa.s:ebola"/><log idref="gammaShape.s:ebola"/>
        {log_clock}
    </logger>
    <logger id="screenlog" spec="Logger" logEvery="{max(L // 50, 10000) if mode == 'mcmc' else 100000}"><log idref="posterior"/><log idref="likelihood"/><log idref="prior"/></logger>
  </mcmc>'''
    head = f'<?xml version="1.0" encoding="UTF-8" standalone="no"?><beast beautitemplate="Standard" namespace="{NS}" version="2.0">\n    <data id="ebola" spec="Alignment" name="alignment">\n{seq}\n    </data>\n'
    maps = '    <map name="Uniform">beast.base.inference.distribution.Uniform</map><map name="Exponential">beast.base.inference.distribution.Exponential</map><map name="LogNormal">beast.base.inference.distribution.LogNormalDistributionModel</map><map name="Gamma">beast.base.inference.distribution.Gamma</map><map name="prior">beast.base.inference.distribution.Prior</map><map name="LaplaceDistribution">beast.base.inference.distribution.LaplaceDistribution</map><map name="Normal">beast.base.inference.distribution.Normal</map>\n'
    if mode == "mcmc": body = mcmc.replace("<mcmc ", "<run ", 1).replace("</mcmc>", "</run>")
    else:
        body = (f'<run spec="modelselection.inference.PathSampler" nrOfSteps="{A.ps_steps}" chainLength="{A.ps_length}" alpha="0.3" rootdir="ps_{kind}" burnInPercentage="20" preBurnin="{A.ps_pre}" deleteOldLogs="true">\n'
                'cd $(dir)\n$(java) -cp $(java.class.path) beast.pkgmgmt.launcher.BeastLauncher $(resume/overwrite) -seed $(seed) beast.xml\n\n' + mcmc + '\n</run>')
    os.makedirs("xml", exist_ok=True); open(out, "w").write(head + maps + body + "\n</beast>\n")
MODELS = ["het_strict", "het_relaxed", "iso", "dr1", "dr2", "dr3"]
if A.cmd == "make":
    for m in MODELS:
        make(m, "mcmc", f"xml/mcmc_{m}.xml"); make(m, "ps", f"xml/ps_{m}.xml")
    print("wrote", len(MODELS) * 2, "XML files; first date", round(first, 3), "last date", round(last, 3))

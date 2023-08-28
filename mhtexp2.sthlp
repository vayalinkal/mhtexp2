{smcl}
{findalias asfradohelp}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "[R] help" "help help"}{...}
{viewerjumpto "Syntax" "mhtexp2##syntax"}{...}
{viewerjumpto "Description" "mhtexp2##description"}{...}
{viewerjumpto "Options" "mhtexp2##options"}{...}
{viewerjumpto "Remarks" "mhtexp2##remarks"}{...}
{title:Title}

{phang}
{bf:mhtexp2} {hline 2} Stata command for the procedure detailed in List, Shaikh, and Vayalinkal (2023). Extends the Stata package mhtexp to allow for covariate-adjustment.


{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:mhtexp2}
{varlist}
{cmd:, } {it:treatment} [{it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opth treatment(varlist)}}treatment status variables {it:varlist}{p_end}
{synopt:{opth subgroup(varname)}}group identifier variable {it:varname}{p_end}
{synopt:{opth controls(varlist)}}covariates to include as controls {it:varlist}{p_end}
{synopt:{opth combo(string)}}compare "treatmentcontrol" or "pairwise"; default is
    {cmd:combo("treatmentcontrol")}{p_end}
{synopt:{opth only(name)}} the numoc*numsub*numpc hypotheses to be tested{p_end}
{synopt:{opth exclude(name)}} the numoc*numsub*numpc hypotheses not to be tested{p_end}
{synopt:{opth bootstrap(integer)}} the number of simulated samples to use{p_end}
{synopt:{opth transitivitycheck(integer)}} whether or not to implement the improvement outlined in Remark 3.8{p_end}
{synoptline}
{p2colreset}{...}

{marker description}{...}
{title:Description}

{pstd}
{cmd:mhtexp2} testing procedure for multiple hypothesis testing that asymptotically controls
familywise error rate and is asymptotically balanced for outcomes specified via {varlist}, while adjusting for baseline covariates{p_end}

{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt treatment(varlist)} user provided variable containing treatment status of the observations, should be integers starting with 0 - if there are many treatments, the command will run much faster with the option transitivitycheck(0) (see below); required.{p_end}

{phang}
{opt subgroup(varname)} user provided variable containing subgroup ids, should be integers starting with 1 (not 0); optional.{p_end}

{phang}
{opt controls(varlist)} user provided list of variables to include as controls; these variables will be incorporated using a "fully interacted" regression specification. 
details provided in List, Shaikh, and Vayalinkal (2021). optional.{p_end}

{phang}
{opt combo(string)} user provided string to specify the comparison between treatment and control.
{cmd:combo("pairwise")} will compare all pairwise comparisons across treatment and control.
The default is {cmd:combo("treatmentcontrol")}, compares each treatment to the control; optional
{p_end}

{phang}
{opt only(name)} N by 3 matrix specifying which hypothesis to be tested; optional.{p_end}

{phang}
{opt exclude(name)} N by 3 matrix specifying which hypothesis not to be tested; optional.{p_end}
{phang}
The matrix in either case should be defined where in each row, column 1 is the outcome, column 2 is the
subgroup and column 3 is the treatment-control comparison. Where...{p_end}
{phang3} 1 <= column 1 <= number of outcomes{p_end}
{phang3} 1 <= column 2 <= number of subgroups{p_end}
{phang3} 1 <= column 3 <= number of treatment-control comparisons{p_end}

{phang}
By default {cmd:mhtexp2} will calculate all hypotheses based on the number of outcomes, subgroups and treatments provided by the user
in {it:varlist} {it:subgroup(varname)} and {it:treatment(varname)}, respectively. Section 4.4 of List, Shaikh and Vayalinkal (2021) simultaniously considers
4 outcome variables, 4 subgroups and 3 treatment conditions, producting a table of 48 hypothesis test. However, there are cases in which you
may only be interested in certain outcome by subgroup by treatment hypothesis. use {opt only} or {opt exclude}.{p_end}


{phang}
{opt bootstrap(integer)} the number of simulated samples. the default is 3000,  but a larger number is recommended when there are a large number of hypotheses; optional.{p_end}


{phang}
{opt transitivitycheck(integer)} whether or not to implement the improvement outlined in Remark 3.8, which can be computationally prohibitive when there are many treatments. Turning this off makes the Remark 3.8 column the same as the Theorem 3.1 column; optional.{p_end}


{marker remarks}{...}
{title:Remarks}

{pstd}
For detailed information on the procedure, see List, Shaikh, and Vayalinkal (2023) (http://individual.utoronto.ca/atom/files/mht2.pdf).{p_end}
{phang2}
{cmd:. ssc install moremata}{p_end}

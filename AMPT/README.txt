AMPT run environment requires AliRoot/latest installed.

TO RUN AMPT:
-> setting number of events:
    goto rungen.C
    modify lines 1 and 2 to number of events desired

-> setting collision energy:
    goto fastMcProduction
    modify line 163 to desired energy in GeV

-> ensure program is AmptDefault
    goto fastMcProduction
    verify line 164 = AmptDefault

-> setting target and projectile
    goto fastMcProduction
    set lines 806, 807 to intended nuclei
    of form ("NUCLEI NAME(if desired)", A, Z)

-> setting Tune
    goto fastMcProduction
    lines 794-833 set tune

-> initiating runs
    load AliPhysics/latest
    $/directory/that/contains/AMPT/programs/> aliroot rungen.C

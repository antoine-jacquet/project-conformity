## Simulations on the evolution of conformity

Consider a group of 2*n* individuals, *n* males and *n* females.

There are 2 diallelic genes. The first one, with alleles *A* and *a*, is expressed in males only and is neutral for fitness. The second one, with alleles *B* and *b*, is expressed in females only: females *b* mate randomly, while females *B* prefer the male phenotype (*A* or *a*) which had higher mating success in the previous generation. For this reason, females *B* are called copier females.

Specifically, call *p* the proportion of matings involving males *A* in the previous generation. Then copier females mate with a male *A* with probability:

*p*^\beta/(*p*^\beta+(1-*p*)^beta),

where \beta is a ‘conformity strength’ parameter (steepness of the conformity curve).

Each female has two children, one male and one female. The recombination rate of the two genes of interest is r (i.e. with probability 1 - r children inherit both alleles from the same parent).

The meta-population is composed of N such groups. In each new generation, a fraction d of the children are selected for migration. The migrants exit their group and are redistributed randomly on the spots left vacant by other migrants (males take male spots and females take female spots to prevent gender bias in groups).


### Initial conditions

Alleles A and a are randomly allocated in the meta-population: in the first generation individuals are attributed either A or a with 50:50 probability. Alleles B and b are also randomly allocated, with an initial proportion q0 of allele B.

### Stopping condition

Simulations stop whenever genetic variation is lost in one of the genes, with a maximum number of generations T.

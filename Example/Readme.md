Example program
====

In this example, we will reproduce figure 1a top from [Mosqueiro &
Maia
PRE](http://neurobiofisica.ifsc.usp.br/Publications/publications.html). Of
course the same program can also reproduce 1a bottom, but this is a
minor change. The idea is to show minimally how to use our functions
and subroutines. To run this example, simply run
```
chmod +x experiment.sh
sh experiment.sh
```

It should take a few minutes to complete. At the end, you should see
two PDF files that should look similar to:

<img src="https://raw.githubusercontent.com/thmosqueiro/Probab-Greenberg-Hastings/master/Example/plot_CDFs.png" width=350px />
<img src="https://raw.githubusercontent.com/thmosqueiro/Probab-Greenberg-Hastings/master/Example/plot_PDFs.png" width=300px />

These results show in green a subcritical avalanche distribution and
in blue, a critical distribution.

First of all, view file size_av.sh as a kind of old-fashioned
Makefile: it will compile the fortran using the environment variables
(see parent folder) and link all libraries. I have removed a bunch of
optimizations that I often use because this might be computer
dependent, and the performance seems already pretty good.

The main file is in sa.f90: you can build almost all experiments with
very similar files. You can also paralelize it with minimal changes in
both sa.f90 and size_sa.sh.

Dependencies
----

To run the example in Example folder, you'll need GnuPlot,
texlife-epstopdf and python 2.7. Pypy is preferable.


Usage examples
----



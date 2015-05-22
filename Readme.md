Probabilistic Greenberg-Hastings
====

In this repository, you'll find an implementation of the probabilistic
Greenberg-Hastings automaton. This implementation is the result of two
years of my PhD with Prof. Dr. Leonardo P. Maia, where we were
investigating neuronal avalanches from the point of view of a simple
toy model that mimics basic features in neurodynamics. This model is
based in some papers ([Haldeman & Beggs -
Phys. Rev.](http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.94.058101#fulltext),
[Kinouchi & Copelli - Nat
Phys.](http://www.nature.com/nphys/journal/v2/n5/full/nphys289.html)). I'm
not just sharing the core code, but also a small example script to
reproduce one of the figures shown in [this
paper](http://journals.aps.org/pre/abstract/10.1103/PhysRevE.88.012712). You
can find this example in Example folder.

Although almost all comments are in portuguese, I plan to eventually
translate it to english someday. This is not a priority though,
especially because anyone can figure things out most of the time. If
you are using this code, please (please!!)  star this repository and
cite [T Mosqueiro & L Maia -- Phys Rev E v88
p012712](http://journals.aps.org/pre/abstract/10.1103/PhysRevE.88.012712).

Dependencies
---

There no real dependencies for the main code (everything that is in
the core folder is self-contained). All you'll need is a fortran
compiler, either Intel's or GCC.

However, to run the example in Example folder, you'll also need
GnuPlot, texlife-epstopdf and python 2.7. Pypy is preferable.


License
---

**tl;dr version:** please, don't sue me and ship a copy of the License
  file with any derived product. Read License file for the actual
  terms.
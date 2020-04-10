**ANANSE**: **AN**\ alysis **A**\ lgorithm for **N**\ etworks **S**\ pecified by **E**\ nhancers
------------------------------------------------------------------------------------------------

What is ANANSE Network?
~~~~~~~~~~~~~~~~~~~~~~~

ANANSE is an analysis framework for key transcription factor prediction
during cell fate switch written in Python. It contains command-line
scripts to predict all TF binding network (``binding``), predict gene
regulatory network (``network``) and infer influence score between two
cell types(tissues) (``influence``). In addition, all this functionality
is available from a Python API.

ANANSE is free and open source research software. If you find it useful
please cite our paper:

-  Xu, Q., Georgiou, G., Veenstra, G.J.C., Zhou, H., and van Heeringen,
   S.J. (2020). **ANANSE: An enhancer network-based computational
   approach for predicting key transcription factors in cell fate
   determination.** https://www.biorxiv.org/

Getting started
~~~~~~~~~~~~~~~

-  The easiest way to install ANANSE is using `bioconda`_ on Linux
   (include `Windows Subsystem for Linux`_ ) or Mac.
-  Have a look at these simple examples to get a taste of what is
   possible.
-  Check out the more detailed tutorials.
-  Full command-line reference can be found here.
-  There’s also an `API documentation`_.

Get help
~~~~~~~~

-  First, check the FAQ for common issues.
-  The preferred way to get support is through the `Github issues
   page`_.
-  Finally, you can reach me by `mail`_ or via `twitter`_.

Full contents
~~~~~~~~~~~~~

-  `Installation`_

   -  `The easiest way to install`_
   -  `Alternative installation`_

-  `Input data`_

   -  `Genome`_
   -  `Motif database`_
   -  `Enhancer data`_
   -  `Expression data`_
   -  `Differential expression data`_

-  `Command-line reference`_

   -  List of tools
   -  Input formats
   -  `Command: ananse binding`_
   -  `Command: ananse network`_
   -  `Command: ananse influence`_

-  `API documentation <API_documentation.md>`__

   -  `Working with ``binding`` Class`_
   -  `Working with ``network`` Class`_
   -  `Working with ``influence`` Class`_

-  `Examples`_

   -  `Command: ananse binding <examples/#motif-database>`__
   -  `Command: ananse network <examples/#motif-database>`__
   -  `Command: ananse influence <examples/#motif-database>`__

-  FAQ

   -  Sorry, motif prediction tool [X] is not supported
   -  I get motifs that have differential scores in

.. _bioconda: https://bioconda.github.io/
.. _Windows Subsystem for Linux: https://docs.microsoft.com/en-us/windows/wsl/install-win10
.. _API documentation: api.md
.. _Github issues page: https://github.com/vanheeringen-lab/ANANSE/issues
.. _mail: mailto:qxuchn@gmail.com
.. _twitter: https://twitter.com/qxuchn
.. _Installation: installation.md
.. _The easiest way to install: installation/#the-easiest-way-to-install
.. _Alternative installation: installation/#alternative-installation
.. _Input data: input_data.md
.. _Genome: input_data/#genome
.. _Motif database: input_data/#motif-database
.. _Enhancer data: input_data/#enhancer-data
.. _Expression data: input_data/#enhancer-data
.. _Differential expression data: input_data/#enhancer-data
.. _Command-line reference: command-line_reference.md
.. _`Command: ananse binding`: command-line_reference/#motif-database
.. _`Command: ananse network`: command-line_reference/#motif-database
.. _`Command: ananse influence`: command-line_reference/#motif-database
.. _Working with ``binding`` Class: API_documentation/#working-with-binding-class
.. _Working with ``network`` Class: API_documentation/#working-with-network-class
.. _Working with ``influence`` Class: API_documentation/#working-with-influence-class
.. _Examples: examples.md
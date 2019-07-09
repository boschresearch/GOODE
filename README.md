# GP BVP solver

This is the companion code for the paper GOODE: A Gaussian Off-the-shelf 
Ordinary Differential Equation Solver by David John, Michael Schober and
Vincent Heuvline. The paper is submitted to ICML 2019 and can be found 
[here](http://proceedings.mlr.press/v97/john19a.html).
The code allows the users to reproduce and extend the results reported in 
the study. Please cite the above paper when reporting, reproducing or 
extending the results.

```
@InProceedings{pmlr-v97-john19a,
  title = 	 {{GOODE}: A {G}aussian Off-The-Shelf Ordinary Differential Equation Solver},
  author = 	 {John, David and Heuveline, Vincent and Schober, Michael},
  booktitle = 	 {Proceedings of the 36th International Conference on Machine Learning},
  pages = 	 {3152--3162},
  year = 	 {2019},
  editor = 	 {Chaudhuri, Kamalika and Salakhutdinov, Ruslan},
  volume = 	 {97},
  series = 	 {Proceedings of Machine Learning Research},
  address = 	 {Long Beach, California, USA},
  month = 	 {09--15 Jun},
  publisher = 	 {PMLR}
}
```

## Purpose of the project

This software is a research prototype, solely developed for and published as
part of the publication cited above. It will neither be maintained nor 
monitored in any way.

## Requirements, how to build, test, install, use, etc.

The code was developed and works with Matlab R2018a. 

Download the bvp Testset from 
https://archimede.dm.uniba.it/~bvpsolvers/testsetbvpsolvers/?page_id=27 
and copy to the folder [external/bvpTestSet](external/bvpTestSet).
Download the code TOM from
https://archimede.dm.uniba.it/~mazzia/mazzia/?page_id=433
and copy to the folder [external/tom](external/tom).

Run the setup.m file to add required folders to the path. The code for the 
GB BVP solver GOODE is in the folder [GPs](GPs). The folder [experiments](experiments) contains
scripts used to produce the results in the paper.   

## License

GP BVP solver is open-sourced under the MIT license. See the
[LICENSE](LICENSE) file for details.

For a list of other open source components included in GP BVP, see the
file [3rd-party-licenses.txt](3rd-party-licenses.txt).


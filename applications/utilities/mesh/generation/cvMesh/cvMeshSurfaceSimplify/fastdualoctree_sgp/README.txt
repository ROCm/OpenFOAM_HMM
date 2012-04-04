
 Statistical optimization of octree searches
 Fast generation of pointerless octree duals
---------------------------------------------------
Rener Castro, Thomas Lewiner, Hélio Lopes, Geovan Tavares, Alex Bordignon
Thomas Lewiner, Vinícius Mello, Adelailson Peixoto, Sinésio Pesco, Hélio Lopes


http://www.matmidia.mat.puc-rio.br/tomlew/publication_page.php?pubkey=octree_cgf
http://www.matmidia.mat.puc-rio.br/tomlew/publication_page.php?pubkey=fastdualoctree_sgp

1. Disclaimer
2. Available code compared to the papers
3. Octrees implementations
4. Morton codes
5. Data
6. Hash table
7. Requirements


Disclaimer
--------------
  This code is implements fast searches and the dual generation of static
or dynamic pointer-based or pointerless octrees. It is provided "as is"
and can be used freely. We just ask the users who would use this code to
mention in their academic paper, patent or software documentation the
original articles:

@article{octree_cgf,
    author = {Rener Castro and Thomas Lewiner and Hélio Lopes and Geovan Tavares and Alex Bordignon},
    title = {Statistical optimization of octree searches},
    year = {2008},
    month = {march},
    journal = {Computer Graphics Forum},
    volume = {27},
    number = {6},
    pages = {1557--1566},
    publisher = {Eurographics},
    doi = {10.1111/j.1467-8659.2007.01104.x},
    url = {\url{http://www.matmidia.mat.puc-rio.br/tomlew/pdfs/octree_cgf.pdf}}
}

@article{fastdualoctree_sgp,
    author = {Thomas Lewiner and Vinícius Mello and Adelailson Peixoto and Sinésio Pesco and Hélio Lopes},
    title = {Fast Generation of Pointerless Octree Duals},
    year = {2010},
    month = {july},
    journal = {Symposium on Geometry Processing 2010},
    booktitle = {Computer Graphics Forum},
    volume = {29},
    number = {5},
    pages = {1661--1669},
    publisher = {Wiley},
    address = {Lyon, France},
    url = {\url{http://www.matmidia.mat.puc-rio.br/tomlew/pdfs/fastdualoctree_sgp.pdf}}
}



If you have any comment, correction, question, please feel free to
contact me at : lewiner@gmail.com .
Thank to David Coeurjolly for the code porting to UNIX.




Available code compared to the paper
----------------------------------------------
  The code contains simple implementations used to generate the result of
our papers. We stripped down the time measurements (which are platform
specific) and the robust dual marching cubes
(http://www.matmidia.mat.puc-rio.br/tomlew/publication_page.php?pubkey=adaptive_implicit_sibgrapi)
code.



Octrees implementations
------------------------------
  The code contains implementations of a pointer-based octree
("ptr_octree.*" files), pointerless octrees with recursive dual
generation ("hash_octree.*"), pointerless octrees with optimized dual
generation for the static strategy representing all the nodes
("opt_octree.*") or only the leaves ("leaf_octree.*"), and the dynamic
strategy ("mem_octree.*"). The octree is chosen at compilation time
through the OCTREE_SWITCH macro.


Morton codes
-----------------
  The Morton code operations in dilated integers are all grouped in the
"morton.*" files.


Data
------
  The data (implicit function or isogrid) is accessed through the
data_access class. The isogrid must be in "iso" format, i.e. a binary
file with 32-bits integers nx, ny and nz for the size of the grid,
32-bits floats xmin, xmax, ymin… for the geometry of the cube, followed
by nx.ny.nz 32-bits floats for the data.


Hash table
-------------
  Finally, we made a very simple naive and simple implementations for the
hashtable, in open ("hash_ptr.h") or close ("has_noptr.h") modes. The
parameters of the hash function are chosen at compilation time through
macros defined in the "hash.h" file.


Requirements
------------

  The graphical interfaces require the openGL libraries with the GLU, glut
and glui extension, available at http://glui.sourceforge.net/ .

  We use Juha Nieminen, Joel Yliluoma's function parser for the implicit
function, and our Marching Cube's implementation
(http://www.matmidia.mat.puc-rio.br/tomlew/publication_page.php?pubkey=marching_cubes_jgt).


Thank you for your interest, and have fun!
http://www.matmidia.mat.puc-rio.br/tomlew/

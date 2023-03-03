  ```
  $ mkdir build
  $ cd build
  $ cmake ..
  $ make
  ```

At this point one can install on this local machine directly as follows:
  ```
  $ sudo make install
  $ sudo ldconfig
  ```
  
Or, generate a .deb file, which can then be transfered to a ppa server so that this library can be installed using the usual apt install sequence.
  ```
  $ sudo cpack -G DEB  <- This will generate a deb file of the version specified in CMakelists.txt in ./_packages/.
  ```
  
One can then install it locally:
   ```
   $ sudo dpkg install ./_packages/*.deb
   ```

Or, better, tranfer it to a ppa repository.  In this case https://github.com/hamilton8415/ppa

To do so, do the following:
  ```
  $ cd ~; git clone git@github.com:hamilton8415/ppa.git
  $ cp ~/FreeSurfaceHydrodynamics/_packages/*.deb ~/ppa/.
  $ cd ~/ppa  
  $ ./scripts/sign_debian.sh
  $ git add [new .deb files]
  $ git commit -a
  $ git push
```

After a few minutes, the new debian file should be availble in the ppa repository, so to install on any machine:
  ```
  $ curl -s --compressed "https://hamilton8415.github.io/ppa/KEY.gpg" | gpg --dearmor | sudo tee /etc/apt/trusted.gpg.d/ppa.gpg >/dev/null
  $ sudo curl -s --compressed -o /etc/apt/sources.list.d/my_list_file.list "https://hamilton8415.github.io/ppa/my_list_file.list"
  $ sudo apt update
  $ sudo apt install libfshydrodynamics 
  ```
  
This should install the relevant files, and run ldconfig so that the library is available




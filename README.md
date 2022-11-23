# tp3-metnum

Antes de intentar compilar el proyecto, recomendamos correr el siguiente script (intentará instalar paquetes en caso de que no se encuentren, por lo que puede llegar a requerir acceso sudo):
$ cd [RUTA_PROYECTO]/src
$ ./install_dependencies.sh
O en otro caso, chequear tener instaladas las siguientes dependencias:
- cmake
- python3
- python notebook
- python pybind11
- python numpy
- python pandas
- python matplotlib
- python scipy
- python networkx
- python seaborn
- python import_ipynb

## INSTRUCCIONES PARA COMPILAR:
$ cd [RUTA_PROYECTO]/src
$ ./compile.sh

(el script de compile se encarga de borrar los archivos innecesarios y de ejecutar el make necesario para compilar)

## INSTRUCCIONES PARA CORRER EL CÓDIGO (luego de compilado):
Todo el código disponible para la cátedra va a estar presente en la notebook de python (gracias a pybind11), por lo que simplemente se debe iniciar la notebook presente en:
$ [RUTA_PROYECTO]/src/Test_tp2.ipynb
Para ejecutar la notebook, es suficiente hacer doble click sobre el archivo (siempre y cuando se tenga instaladas las dependencias mencionadas anteriormente)

## INSTRUCCIONES PARA LIMPIAR LOS ARCHIVOS CREADOS AL COMPILAR:
$ cd [RUTA_PROYECTO]/src
$ ./clear_build.sh
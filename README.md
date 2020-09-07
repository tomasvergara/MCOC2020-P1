# MCOC2020-P1
## Entrega 1
![grafico_balistica](https://user-images.githubusercontent.com/69252038/91106071-c7f0da00-e63f-11ea-9bcd-77b35d22a18f.png)

* Grafico de balistica con diferentes valores de velocidad del viento
## Entrega 5
### Parte 1
![grafico_xyz](https://user-images.githubusercontent.com/69252038/92355789-5b53f180-f0bb-11ea-9d30-6e9d4b9706bd.png)
![grafico_deriva](https://user-images.githubusercontent.com/69252038/92371914-9e6d8f00-f0d2-11ea-981f-f4e2fe9086f9.png)
* En el grafico de la izquierda se observan los comportamientos que tiene el satelite con la orbita predicha "basica" y la real entregada en el archivo EOF. En el grafico de la derecha se observa que la diferencia final entre ambas orbitas calculadas es de: 597,9 km. 
### Parte 2
![ode_vs_eul_deriva_nsub=1](https://user-images.githubusercontent.com/69252038/92357393-21381f00-f0be-11ea-93fc-31a995ebedfa.png)
* En el grafico se observa la deriva con un numero de subdivisiones igual a 1, con las distintas soluciones de "odeint" y "eulerint". Se ve que "eulerint"(color azul) demuestra un comportamiento oscilante que se aleja cada vez mas de la Tierra, siendo muy poco efectivo en su desempeño, tanto asi, que la deriva entre ambas soluciones al momento final del tiempo es de 19.689,96 km. 
* Los tiempos de realizacion de las soluciones son: Odeint = 0.30618 segundos , Eulerint = 0.99744 segundos
### Parte 3
![ode_vs_eul_deriva_nsub=1700_parte2y3](https://user-images.githubusercontent.com/69252038/92373845-4e43fc00-f0d5-11ea-85eb-29d3ed15d202.png)
* Al realizar varias iteraciones para analizar el comportamiento de la funcion "Eulerint", se vio que al sobrepasar el Nsub = 1000 se veia un porcentaje de error cercano al buscado. Hasta que se realizó con un Nsub = 1700, el cual arrojó el resultado mas positivo de los experimentados (1000, 1200 y 1500). La iteracion con Nsub = 1700, con Eulerint tardó 2022,356 segundos, los cuales son cerca de 33 min de demora. Al momento que python mostró el grafico, vi que en el codigo habia establecido ciertos limites que no permiten visualizar bien el comportamiento de ambas curvas, eso si, se puede ver que la linea azul esta muy cerca del 0, esto significa que el error es muy pequeño. 
### Parte 4
#### Con J2:
![grafico_xyz_J2](https://user-images.githubusercontent.com/69252038/92376472-045d1500-f0d9-11ea-9a03-44f6f9e4c731.png)
![grafico_deriva_J2](https://user-images.githubusercontent.com/69252038/92376492-0b842300-f0d9-11ea-8362-0c668803a13d.png)
* Se puede observar que el grafico de la izquierda no hay grandes cambios visuales con respecto a los obtenidos en la Parte 1, pero con el grafico de la derecha queda demostrado que la incorporacion de J2 en la funcion ayudó a disminuir en un 50% nuestro error aproximadamente con un Delta = 319.89 km, lo cual es bueno, pero no lo suficiente como ocurrió con los graficos presentados en Canvas por el profesor.
![ode_vs_eul_deriva_nsub=1_J2](https://user-images.githubusercontent.com/69252038/92376521-12129a80-f0d9-11ea-895a-2fe583e17266.png)
* En este grafico tambien queda demostrado que hay mejores en el rendimiento del codigo, ya que la curva roja (Odeint) disminuye su "pendiente" y eso indica que el calculo mediante esta solucion esta mas cerca que sin J2.
* Este codigo se demoro en correr 3.793 segundos
#### Con J2 y J3:
![grafico_xyz_J2yJ3](https://user-images.githubusercontent.com/69252038/92377686-e5f81900-f0da-11ea-8080-0b95ff176d80.png)
![grafico_deriva_J2yJ3](https://user-images.githubusercontent.com/69252038/92377698-eb556380-f0da-11ea-83ba-899b6980e8c2.png)
* Se puede observar que el grafico de la izquierda no hay grandes cambios visuales con respecto a los obtenidos en la Parte 1, pero con el grafico de la derecha se ve que la incorporacion de J3 ayudó a disminuir de manera leve el error con solo J2 con un Delta = 317.59 km, lo cual es bueno, ya que se respeta lo presentado en los graficos de Canvas. 
![ode_vs_eul_deriva_nsub=1_J2yJ3](https://user-images.githubusercontent.com/69252038/92377712-f4decb80-f0da-11ea-941e-dfef0dab5bd6.png)
* En este grafico, no se nota visualmente una diferencia con el anterior, pero si hay un diferencias en los calculos, viendose reflejado en un disminucion en el punto mas alto de la curva roja, lo que minimiza el error.  
* Curiosamente, este codigo se demoro en correr 3.657 segundos, siendo menor al tiempo empleado con solo J2 incluido. 

# rollercoaster-pathfinder
## This is a prelimiary version of the code

This code reconstructs a full rollercoaster track from the linear and angular accelerations measured during the ride. The data was taken using the Roller Coaster mode in the Physics Toolbox Suite mobile app. This data is stored in the .csv file. First, all points are moved to the same system of reference, based at the initial point. Calculations are made based on a self-developed integrating method. With this, the relative angle at each time is calculated from the angular velocity. Linear velocity and position are similarly calculated from linear acceleration. A coordinate redefinition is made since the position of the mobile phone (the data recorder) changes during the ride. As a check, the initial and final positions should be the same.

---

Aquest codi pretén reconstruir el camí d'una muntanya russa a partir de l'acceleració lineal i la velocitat angular mesurades durant el recorregut. Les dades s'han pres amb el mode Roller Coaster de l'aplicació mòbil Physics Toolbox Suite. Les dades preses es troben al fitxer .csv. Al codi, primerament tots els punts es mouen al mateix sistema de referència, centrar al punt inicial. Els càlculs es basen en un mètode d'integració propi, adaptat a les dades que es tenen. Amb aquest mètode es troba l'angle relatiu a cada punt a partir de la velocitat angular. Similarment es calculen la velocitat lineal i la posició a partir de l'acceleració angular. Com a comprovació, les posicions inicial i final haurien de ser les mateixes.

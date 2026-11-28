clear; close; clc;
device = ThermodynamicDevice("diffuser")
state0 = ThermodynamicState("freeStream",["ideal gas","caloricly perfect"],"previousDevice",device)
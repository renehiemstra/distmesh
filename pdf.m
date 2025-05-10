clear; close all; clc;

x = [-3:.1:3];
y = @(x) normpdf(x,0,1) ./  normpdf(0,0,1);

plot(x,y(x))

y(0)
y(1)
y(2)
y(3)
y(4)
y(5)
y(6)

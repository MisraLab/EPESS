N=100;
points_1 = zeros(length(t1),2);
points_2 = zeros(length(t2),2);
points_3 = zeros(length(angle_slice),2);
points_4 = zeros(length(angle_slice1),2);

prop_0 = zeros(N,2);
prop_1 = zeros(N,2);
prop_2 = zeros(N,2);
prop_3 = zeros(N,2);
prop_4 = zeros(N,2);
prop_5 = zeros(N,2);
prop_6 = zeros(N,2);

theta_0 = 0:2*pi/N:2*pi;
theta_1 = angle_slice(1):angle_slice(2)/N:angle_slice(2);
theta_2 = angle_slice(3):(angle_slice(4)-angle_slice(3))/N:angle_slice(4);
theta_3 = angle_slice(5):(angle_slice(6)-angle_slice(5))/N:angle_slice(6);


theta_4 = -2*pi:(angle_slice1(1)+2*pi)/N:angle_slice1(1);
theta_5 = angle_slice1(2):(angle_slice1(3)-angle_slice1(2))/N:angle_slice1(3);
theta_6 = angle_slice1(4):(0 - angle_slice1(4))/N:0;



lB=[10,-1];
uB=[11, 1];

xx_1 = b;
nu_1=a;

for i=1:length(t1)
    points_1(i, :) = xx_1*cos(t1(i)) + nu_1*sin(t1(i));
    points_2(i, :) = xx_1*cos(t2(i)) + nu_1*sin(t2(i));

end

for i=1:length(angle_slice)
    points_3(i, :) = xx_1*cos(angle_slice(i)) + nu_1*sin(angle_slice(i));
    points_4(i, :) = xx_1*cos(angle_slice1(i)) + nu_1*sin(angle_slice1(i));

end

for i =1:N
    prop_0(i,:) = xx_1*cos(theta_0(i)) + nu_1*sin(theta_0(i));
    prop_1(i,:) = xx_1*cos(theta_1(i)) + nu_1*sin(theta_1(i));
    prop_2(i,:) = xx_1*cos(theta_2(i)) + nu_1*sin(theta_2(i));
    prop_3(i,:) = xx_1*cos(theta_3(i)) + nu_1*sin(theta_3(i));
    prop_4(i,:) = xx_1*cos(theta_4(i)) + nu_1*sin(theta_4(i));
    prop_5(i,:) = xx_1*cos(theta_5(i)) + nu_1*sin(theta_5(i));
    prop_6(i,:) = xx_1*cos(theta_6(i)) + nu_1*sin(theta_6(i));
   
end

line(prop_0(:,1), prop_0(:,2))
hold on 
line(prop_1(:,1), prop_1(:,2))
hold on 
line(prop_2(:,1), prop_2(:,2))
hold on
line(prop_3(:,1), prop_3(:,2))
hold on
line(prop_4(:,1), prop_4(:,2))
hold on
line(prop_5(:,1), prop_5(:,2))
hold on
line(prop_6(:,1), prop_6(:,2))
hold on
plot(xx_1(1),xx_1(2),'o')
hold on
plot(nu_1(1), nu_1(2), 'o')
hold on
plot(xx_1(1), xx_1(2), '*')
hold on
plot(xx_2(1), xx_2(2), '*')
hold on
plot(xx_3(1), xx_3(2), '*')
hold on
plot(xx_4(1), xx_4(2), '*')
hold on
rectangle('Position',[10 -1 1 2])
hold off
axis([lB(1)-2, uB(1)+2, lB(2)-1, uB(2)+1])
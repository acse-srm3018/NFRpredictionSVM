function dy = ODE20(t,y)
global C20 A20 la om uu20
u=uu20;
C=C20;
A=A20;
dy = zeros(2*(length(u)-1),1);% a column vector
dy(1) = (1/om)*(C(1,1)*y(1)+C(1,2)*y(2)+C(1,3)*y(3)+C(1,4)*y(4)+C(1,5)*y(5)+C(1,6)*y(6)+C(1,7)*y(7)+C(1,8)*y(8)+C(1,9)*y(9)+C(1,10)*y(10)...
  +C(1,11)*y(11)+C(1,12)*y(12)+C(1,13)*y(13)+C(1,14)*y(14)+C(1,15)*y(15)+C(1,16)*y(16)+C(1,17)*y(17)+C(1,18)*y(18)+C(1,19)*y(19)+C(1,20)*y(20)-la*(y(1)-y(21))...
    +(C(1,end)/A(end,end))*((-0.5/(u(end)-(u(end)^0.5)))-(A(end,1)*y(1)+A(end,2)*y(2)+A(end,3)*y(3)+A(end,4)*y(4)+A(end,5)*y(5)+A(end,6)*y(6)+A(end,7)*y(7)+A(end,8)*y(8)+A(end,9)*y(9)+A(end,10)*y(10)...
    +A(end,11)*y(11)+A(end,12)*y(12)+A(end,13)*y(13)+A(end,14)*y(14)+A(end,15)*y(15)+A(end,16)*y(16)+A(end,17)*y(17)+A(end,18)*y(18)+A(end,19)*y(19)+A(end,20)*y(20))));

dy(2) = (1/om)*(C(2,1)*y(1)+C(2,2)*y(2)+C(2,3)*y(3)+C(2,4)*y(4)+C(2,5)*y(5)+C(2,6)*y(6)+C(2,7)*y(7)+C(2,8)*y(8)+C(2,9)*y(9)+C(2,10)*y(10)...
  +C(2,11)*y(11)+C(2,12)*y(12)+C(2,13)*y(13)+C(2,14)*y(14)+C(2,15)*y(15)+C(2,16)*y(16)+C(2,17)*y(17)+C(2,18)*y(18)+C(2,19)*y(19)+C(2,20)*y(20)-la*(y(2)-y(22))...
    +(C(2,end)/A(end,end))*((-0.5/(u(end)-(u(end)^0.5)))-(A(end,1)*y(1)+A(end,2)*y(2)+A(end,3)*y(3)+A(end,4)*y(4)+A(end,5)*y(5)+A(end,6)*y(6)+A(end,7)*y(7)+A(end,8)*y(8)+A(end,9)*y(9)+A(end,10)*y(10)...
    +A(end,11)*y(11)+A(end,12)*y(12)+A(end,13)*y(13)+A(end,14)*y(14)+A(end,15)*y(15)+A(end,16)*y(16)+A(end,17)*y(17)+A(end,18)*y(18)+A(end,19)*y(19)+A(end,20)*y(20))));


dy(3) = (1/om)*(C(3,1)*y(1)+C(3,2)*y(2)+C(3,3)*y(3)+C(3,4)*y(4)+C(3,5)*y(5)+C(3,6)*y(6)+C(3,7)*y(7)+C(3,8)*y(8)+C(3,9)*y(9)+C(3,10)*y(10)...
  +C(3,11)*y(11)+C(3,12)*y(12)+C(3,13)*y(13)+C(3,14)*y(14)+C(3,15)*y(15)+C(3,16)*y(16)+C(3,17)*y(17)+C(3,18)*y(18)+C(3,19)*y(19)+C(3,20)*y(20)-la*(y(3)-y(23))...
    +(C(3,end)/A(end,end))*((-0.5/(u(end)-(u(end)^0.5)))-(A(end,1)*y(1)+A(end,2)*y(2)+A(end,3)*y(3)+A(end,4)*y(4)+A(end,5)*y(5)+A(end,6)*y(6)+A(end,7)*y(7)+A(end,8)*y(8)+A(end,9)*y(9)+A(end,10)*y(10)...
    +A(end,11)*y(11)+A(end,12)*y(12)+A(end,13)*y(13)+A(end,14)*y(14)+A(end,15)*y(15)+A(end,16)*y(16)+A(end,17)*y(17)+A(end,18)*y(18)+A(end,19)*y(19)+A(end,20)*y(20))));

dy(4)=(1/om)*(C(4,1)*y(1)+C(4,2)*y(2)+C(4,3)*y(3)+C(4,4)*y(4)+C(4,5)*y(5)+C(4,6)*y(6)+C(4,7)*y(7)+C(4,8)*y(8)+C(4,9)*y(9)+C(4,10)*y(10)...
  +C(4,11)*y(11)+C(4,12)*y(12)+C(4,13)*y(13)+C(4,14)*y(14)+C(4,15)*y(15)+C(4,16)*y(16)+C(4,17)*y(17)+C(4,18)*y(18)+C(4,19)*y(19)+C(4,20)*y(20)-la*(y(4)-y(24))...
    +(C(4,end)/A(end,end))*((-0.5/(u(end)-(u(end)^0.5)))-(A(end,1)*y(1)+A(end,2)*y(2)+A(end,3)*y(3)+A(end,4)*y(4)+A(end,5)*y(5)+A(end,6)*y(6)+A(end,7)*y(7)+A(end,8)*y(8)+A(end,9)*y(9)+A(end,10)*y(10)...
    +A(end,11)*y(11)+A(end,12)*y(12)+A(end,13)*y(13)+A(end,14)*y(14)+A(end,15)*y(15)+A(end,16)*y(16)+A(end,17)*y(17)+A(end,18)*y(18)+A(end,19)*y(19)+A(end,20)*y(20))));

dy(5)=(1/om)*(C(5,1)*y(1)+C(5,2)*y(2)+C(5,3)*y(3)+C(5,4)*y(4)+C(5,5)*y(5)+C(5,6)*y(6)+C(5,7)*y(7)+C(5,8)*y(8)+C(5,9)*y(9)+C(5,10)*y(10)...
  +C(5,11)*y(11)+C(5,12)*y(12)+C(5,13)*y(13)+C(5,14)*y(14)+C(5,15)*y(15)+C(5,16)*y(16)+C(5,17)*y(17)+C(5,18)*y(18)+C(5,19)*y(19)+C(5,20)*y(20)-la*(y(5)-y(25))...
    +(C(5,end)/A(end,end))*((-0.5/(u(end)-(u(end)^0.5)))-(A(end,1)*y(1)+A(end,2)*y(2)+A(end,3)*y(3)+A(end,4)*y(4)+A(end,5)*y(5)+A(end,6)*y(6)+A(end,7)*y(7)+A(end,8)*y(8)+A(end,9)*y(9)+A(end,10)*y(10)...
    +A(end,11)*y(11)+A(end,12)*y(12)+A(end,13)*y(13)+A(end,14)*y(14)+A(end,15)*y(15)+A(end,16)*y(16)+A(end,17)*y(17)+A(end,18)*y(18)+A(end,19)*y(19)+A(end,20)*y(20))));

dy(6)=(1/om)*(C(6,1)*y(1)+C(6,2)*y(2)+C(6,3)*y(3)+C(6,4)*y(4)+C(6,5)*y(5)+C(6,6)*y(6)+C(6,7)*y(7)+C(6,8)*y(8)+C(6,9)*y(9)+C(6,10)*y(10)...
  +C(6,11)*y(11)+C(6,12)*y(12)+C(6,13)*y(13)+C(6,14)*y(14)+C(6,15)*y(15)+C(6,16)*y(16)+C(6,17)*y(17)+C(6,18)*y(18)+C(6,19)*y(19)+C(6,20)*y(20)-la*(y(6)-y(26))...
    +(C(6,end)/A(end,end))*((-0.5/(u(end)-(u(end)^0.5)))-(A(end,1)*y(1)+A(end,2)*y(2)+A(end,3)*y(3)+A(end,4)*y(4)+A(end,5)*y(5)+A(end,6)*y(6)+A(end,7)*y(7)+A(end,8)*y(8)+A(end,9)*y(9)+A(end,10)*y(10)...
    +A(end,11)*y(11)+A(end,12)*y(12)+A(end,13)*y(13)+A(end,14)*y(14)+A(end,15)*y(15)+A(end,16)*y(16)+A(end,17)*y(17)+A(end,18)*y(18)+A(end,19)*y(19)+A(end,20)*y(20))));

dy(7)=(1/om)*(C(7,1)*y(1)+C(7,2)*y(2)+C(7,3)*y(3)+C(7,4)*y(4)+C(7,5)*y(5)+C(7,6)*y(6)+C(7,7)*y(7)+C(7,8)*y(8)+C(7,9)*y(9)+C(7,10)*y(10)...
  +C(7,11)*y(11)+C(7,12)*y(12)+C(7,13)*y(13)+C(7,14)*y(14)+C(7,15)*y(15)+C(7,16)*y(16)+C(7,17)*y(17)+C(7,18)*y(18)+C(7,19)*y(19)+C(7,20)*y(20)-la*(y(7)-y(27))...
    +(C(7,end)/A(end,end))*((-0.5/(u(end)-(u(end)^0.5)))-(A(end,1)*y(1)+A(end,2)*y(2)+A(end,3)*y(3)+A(end,4)*y(4)+A(end,5)*y(5)+A(end,6)*y(6)+A(end,7)*y(7)+A(end,8)*y(8)+A(end,9)*y(9)+A(end,10)*y(10)...
    +A(end,11)*y(11)+A(end,12)*y(12)+A(end,13)*y(13)+A(end,14)*y(14)+A(end,15)*y(15)+A(end,16)*y(16)+A(end,17)*y(17)+A(end,18)*y(18)+A(end,19)*y(19)+A(end,20)*y(20))));

dy(8)=(1/om)*(C(8,1)*y(1)+C(8,2)*y(2)+C(8,3)*y(3)+C(8,4)*y(4)+C(8,5)*y(5)+C(8,6)*y(6)+C(8,7)*y(7)+C(8,8)*y(8)+C(8,9)*y(9)+C(8,10)*y(10)...
  +C(8,11)*y(11)+C(8,12)*y(12)+C(8,13)*y(13)+C(8,14)*y(14)+C(8,15)*y(15)+C(8,16)*y(16)+C(8,17)*y(17)+C(8,18)*y(18)+C(8,19)*y(19)+C(8,20)*y(20)-la*(y(8)-y(28))...
    +(C(8,end)/A(end,end))*((-0.5/(u(end)-(u(end)^0.5)))-(A(end,1)*y(1)+A(end,2)*y(2)+A(end,3)*y(3)+A(end,4)*y(4)+A(end,5)*y(5)+A(end,6)*y(6)+A(end,7)*y(7)+A(end,8)*y(8)+A(end,9)*y(9)+A(end,10)*y(10)...
    +A(end,11)*y(11)+A(end,12)*y(12)+A(end,13)*y(13)+A(end,14)*y(14)+A(end,15)*y(15)+A(end,16)*y(16)+A(end,17)*y(17)+A(end,18)*y(18)+A(end,19)*y(19)+A(end,20)*y(20))));

dy(9)=(1/om)*(C(9,1)*y(1)+C(9,2)*y(2)+C(9,3)*y(3)+C(9,4)*y(4)+C(9,5)*y(5)+C(9,6)*y(6)+C(9,7)*y(7)+C(9,8)*y(8)+C(9,9)*y(9)+C(9,10)*y(10)...
  +C(9,11)*y(11)+C(9,12)*y(12)+C(9,13)*y(13)+C(9,14)*y(14)+C(9,15)*y(15)+C(9,16)*y(16)+C(9,17)*y(17)+C(9,18)*y(18)+C(9,19)*y(19)+C(9,20)*y(20)-la*(y(9)-y(29))...
    +(C(9,end)/A(end,end))*((-0.5/(u(end)-(u(end)^0.5)))-(A(end,1)*y(1)+A(end,2)*y(2)+A(end,3)*y(3)+A(end,4)*y(4)+A(end,5)*y(5)+A(end,6)*y(6)+A(end,7)*y(7)+A(end,8)*y(8)+A(end,9)*y(9)+A(end,10)*y(10)...
    +A(end,11)*y(11)+A(end,12)*y(12)+A(end,13)*y(13)+A(end,14)*y(14)+A(end,15)*y(15)+A(end,16)*y(16)+A(end,17)*y(17)+A(end,18)*y(18)+A(end,19)*y(19)+A(end,20)*y(20))));

dy(10)=(1/om)*(C(10,1)*y(1)+C(10,2)*y(2)+C(10,3)*y(3)+C(10,4)*y(4)+C(10,5)*y(5)+C(10,6)*y(6)+C(10,7)*y(7)+C(10,8)*y(8)+C(10,9)*y(9)+C(10,10)*y(10)...
  +C(10,11)*y(11)+C(10,12)*y(12)+C(10,13)*y(13)+C(10,14)*y(14)+C(10,15)*y(15)+C(10,16)*y(16)+C(10,17)*y(17)+C(10,18)*y(18)+C(10,19)*y(19)+C(10,20)*y(20)-la*(y(10)-y(30))...
    +(C(10,end)/A(end,end))*((-0.5/(u(end)-(u(end)^0.5)))-(A(end,1)*y(1)+A(end,2)*y(2)+A(end,3)*y(3)+A(end,4)*y(4)+A(end,5)*y(5)+A(end,6)*y(6)+A(end,7)*y(7)+A(end,8)*y(8)+A(end,9)*y(9)+A(end,10)*y(10)...
    +A(end,11)*y(11)+A(end,12)*y(12)+A(end,13)*y(13)+A(end,14)*y(14)+A(end,15)*y(15)+A(end,16)*y(16)+A(end,17)*y(17)+A(end,18)*y(18)+A(end,19)*y(19)+A(end,20)*y(20))));

dy(11)=(1/om)*(C(11,1)*y(1)+C(11,2)*y(2)+C(11,3)*y(3)+C(11,4)*y(4)+C(11,5)*y(5)+C(11,6)*y(6)+C(11,7)*y(7)+C(11,8)*y(8)+C(11,9)*y(9)+C(11,10)*y(10)...
  +C(11,11)*y(11)+C(11,12)*y(12)+C(11,13)*y(13)+C(11,14)*y(14)+C(11,15)*y(15)+C(11,16)*y(16)+C(11,17)*y(17)+C(11,18)*y(18)+C(11,19)*y(19)+C(11,20)*y(20)-la*(y(11)-y(31))...
    +(C(11,end)/A(end,end))*((-0.5/(u(end)-(u(end)^0.5)))-(A(end,1)*y(1)+A(end,2)*y(2)+A(end,3)*y(3)+A(end,4)*y(4)+A(end,5)*y(5)+A(end,6)*y(6)+A(end,7)*y(7)+A(end,8)*y(8)+A(end,9)*y(9)+A(end,10)*y(10)...
    +A(end,11)*y(11)+A(end,12)*y(12)+A(end,13)*y(13)+A(end,14)*y(14)+A(end,15)*y(15)+A(end,16)*y(16)+A(end,17)*y(17)+A(end,18)*y(18)+A(end,19)*y(19)+A(end,20)*y(20))));

dy(12)=(1/om)*(C(12,1)*y(1)+C(12,2)*y(2)+C(12,3)*y(3)+C(12,4)*y(4)+C(12,5)*y(5)+C(12,6)*y(6)+C(12,7)*y(7)+C(12,8)*y(8)+C(12,9)*y(9)+C(12,10)*y(10)...
  +C(12,11)*y(11)+C(12,12)*y(12)+C(12,13)*y(13)+C(12,14)*y(14)+C(12,15)*y(15)+C(12,16)*y(16)+C(12,17)*y(17)+C(12,18)*y(18)+C(12,19)*y(19)+C(12,20)*y(20)-la*(y(12)-y(32))...
    +(C(12,end)/A(end,end))*((-0.5/(u(end)-(u(end)^0.5)))-(A(end,1)*y(1)+A(end,2)*y(2)+A(end,3)*y(3)+A(end,4)*y(4)+A(end,5)*y(5)+A(end,6)*y(6)+A(end,7)*y(7)+A(end,8)*y(8)+A(end,9)*y(9)+A(end,10)*y(10)...
    +A(end,11)*y(11)+A(end,12)*y(12)+A(end,13)*y(13)+A(end,14)*y(14)+A(end,15)*y(15)+A(end,16)*y(16)+A(end,17)*y(17)+A(end,18)*y(18)+A(end,19)*y(19)+A(end,20)*y(20))));

dy(13)=(1/om)*(C(13,1)*y(1)+C(13,2)*y(2)+C(13,3)*y(3)+C(13,4)*y(4)+C(13,5)*y(5)+C(13,6)*y(6)+C(13,7)*y(7)+C(13,8)*y(8)+C(13,9)*y(9)+C(13,10)*y(10)...
  +C(13,11)*y(11)+C(13,12)*y(12)+C(13,13)*y(13)+C(13,14)*y(14)+C(13,15)*y(15)+C(13,16)*y(16)+C(13,17)*y(17)+C(13,18)*y(18)+C(13,19)*y(19)+C(13,20)*y(20)-la*(y(13)-y(33))...
    +(C(13,end)/A(end,end))*((-0.5/(u(end)-(u(end)^0.5)))-(A(end,1)*y(1)+A(end,2)*y(2)+A(end,3)*y(3)+A(end,4)*y(4)+A(end,5)*y(5)+A(end,6)*y(6)+A(end,7)*y(7)+A(end,8)*y(8)+A(end,9)*y(9)+A(end,10)*y(10)...
    +A(end,11)*y(11)+A(end,12)*y(12)+A(end,13)*y(13)+A(end,14)*y(14)+A(end,15)*y(15)+A(end,16)*y(16)+A(end,17)*y(17)+A(end,18)*y(18)+A(end,19)*y(19)+A(end,20)*y(20))));

dy(14)=(1/om)*(C(14,1)*y(1)+C(14,2)*y(2)+C(14,3)*y(3)+C(14,4)*y(4)+C(14,5)*y(5)+C(14,6)*y(6)+C(14,7)*y(7)+C(14,8)*y(8)+C(14,9)*y(9)+C(14,10)*y(10)...
  +C(14,11)*y(11)+C(14,12)*y(12)+C(14,13)*y(13)+C(14,14)*y(14)+C(14,15)*y(15)+C(14,16)*y(16)+C(14,17)*y(17)+C(14,18)*y(18)+C(14,19)*y(19)+C(14,20)*y(20)-la*(y(14)-y(34))...
    +(C(14,end)/A(end,end))*((-0.5/(u(end)-(u(end)^0.5)))-(A(end,1)*y(1)+A(end,2)*y(2)+A(end,3)*y(3)+A(end,4)*y(4)+A(end,5)*y(5)+A(end,6)*y(6)+A(end,7)*y(7)+A(end,8)*y(8)+A(end,9)*y(9)+A(end,10)*y(10)...
    +A(end,11)*y(11)+A(end,12)*y(12)+A(end,13)*y(13)+A(end,14)*y(14)+A(end,15)*y(15)+A(end,16)*y(16)+A(end,17)*y(17)+A(end,18)*y(18)+A(end,19)*y(19)+A(end,20)*y(20))));

dy(15)=(1/om)*(C(15,1)*y(1)+C(15,2)*y(2)+C(15,3)*y(3)+C(15,4)*y(4)+C(15,5)*y(5)+C(15,6)*y(6)+C(15,7)*y(7)+C(15,8)*y(8)+C(15,9)*y(9)+C(15,10)*y(10)...
  +C(15,11)*y(11)+C(15,12)*y(12)+C(15,13)*y(13)+C(15,14)*y(14)+C(15,15)*y(15)+C(15,16)*y(16)+C(15,17)*y(17)+C(15,18)*y(18)+C(15,19)*y(19)+C(15,20)*y(20)-la*(y(15)-y(35))...
    +(C(15,end)/A(end,end))*((-0.5/(u(end)-(u(end)^0.5)))-(A(end,1)*y(1)+A(end,2)*y(2)+A(end,3)*y(3)+A(end,4)*y(4)+A(end,5)*y(5)+A(end,6)*y(6)+A(end,7)*y(7)+A(end,8)*y(8)+A(end,9)*y(9)+A(end,10)*y(10)...
    +A(end,11)*y(11)+A(end,12)*y(12)+A(end,13)*y(13)+A(end,14)*y(14)+A(end,15)*y(15)+A(end,16)*y(16)+A(end,17)*y(17)+A(end,18)*y(18)+A(end,19)*y(19)+A(end,20)*y(20))));

dy(16)=(1/om)*(C(16,1)*y(1)+C(16,2)*y(2)+C(16,3)*y(3)+C(16,4)*y(4)+C(16,5)*y(5)+C(16,6)*y(6)+C(16,7)*y(7)+C(16,8)*y(8)+C(16,9)*y(9)+C(16,10)*y(10)...
  +C(16,11)*y(11)+C(16,12)*y(12)+C(16,13)*y(13)+C(16,14)*y(14)+C(16,15)*y(15)+C(16,16)*y(16)+C(16,17)*y(17)+C(16,18)*y(18)+C(16,19)*y(19)+C(16,20)*y(20)-la*(y(16)-y(36))...
    +(C(16,end)/A(end,end))*((-0.5/(u(end)-(u(end)^0.5)))-(A(end,1)*y(1)+A(end,2)*y(2)+A(end,3)*y(3)+A(end,4)*y(4)+A(end,5)*y(5)+A(end,6)*y(6)+A(end,7)*y(7)+A(end,8)*y(8)+A(end,9)*y(9)+A(end,10)*y(10)...
    +A(end,11)*y(11)+A(end,12)*y(12)+A(end,13)*y(13)+A(end,14)*y(14)+A(end,15)*y(15)+A(end,16)*y(16)+A(end,17)*y(17)+A(end,18)*y(18)+A(end,19)*y(19)+A(end,20)*y(20))));

dy(17)=(1/om)*(C(17,1)*y(1)+C(17,2)*y(2)+C(17,3)*y(3)+C(17,4)*y(4)+C(17,5)*y(5)+C(17,6)*y(6)+C(17,7)*y(7)+C(17,8)*y(8)+C(17,9)*y(9)+C(17,10)*y(10)...
  +C(17,11)*y(11)+C(17,12)*y(12)+C(17,13)*y(13)+C(17,14)*y(14)+C(17,15)*y(15)+C(17,16)*y(16)+C(17,17)*y(17)+C(17,18)*y(18)+C(17,19)*y(19)+C(17,20)*y(20)-la*(y(17)-y(37))...
    +(C(17,end)/A(end,end))*((-0.5/(u(end)-(u(end)^0.5)))-(A(end,1)*y(1)+A(end,2)*y(2)+A(end,3)*y(3)+A(end,4)*y(4)+A(end,5)*y(5)+A(end,6)*y(6)+A(end,7)*y(7)+A(end,8)*y(8)+A(end,9)*y(9)+A(end,10)*y(10)...
    +A(end,11)*y(11)+A(end,12)*y(12)+A(end,13)*y(13)+A(end,14)*y(14)+A(end,15)*y(15)+A(end,16)*y(16)+A(end,17)*y(17)+A(end,18)*y(18)+A(end,19)*y(19)+A(end,20)*y(20))));

dy(18)=(1/om)*(C(18,1)*y(1)+C(18,2)*y(2)+C(18,3)*y(3)+C(18,4)*y(4)+C(18,5)*y(5)+C(18,6)*y(6)+C(18,7)*y(7)+C(18,8)*y(8)+C(18,9)*y(9)+C(18,10)*y(10)...
  +C(18,11)*y(11)+C(18,12)*y(12)+C(18,13)*y(13)+C(18,14)*y(14)+C(18,15)*y(15)+C(18,16)*y(16)+C(18,17)*y(17)+C(18,18)*y(18)+C(18,19)*y(19)+C(18,20)*y(20)-la*(y(18)-y(38))...
    +(C(18,end)/A(end,end))*((-0.5/(u(end)-(u(end)^0.5)))-(A(end,1)*y(1)+A(end,2)*y(2)+A(end,3)*y(3)+A(end,4)*y(4)+A(end,5)*y(5)+A(end,6)*y(6)+A(end,7)*y(7)+A(end,8)*y(8)+A(end,9)*y(9)+A(end,10)*y(10)...
    +A(end,11)*y(11)+A(end,12)*y(12)+A(end,13)*y(13)+A(end,14)*y(14)+A(end,15)*y(15)+A(end,16)*y(16)+A(end,17)*y(17)+A(end,18)*y(18)+A(end,19)*y(19)+A(end,20)*y(20))));

dy(19)=(1/om)*(C(19,1)*y(1)+C(19,2)*y(2)+C(19,3)*y(3)+C(19,4)*y(4)+C(19,5)*y(5)+C(19,6)*y(6)+C(19,7)*y(7)+C(19,8)*y(8)+C(19,9)*y(9)+C(19,10)*y(10)...
  +C(19,11)*y(11)+C(19,12)*y(12)+C(19,13)*y(13)+C(19,14)*y(14)+C(19,15)*y(15)+C(19,16)*y(16)+C(19,17)*y(17)+C(19,18)*y(18)+C(19,19)*y(19)+C(19,20)*y(20)-la*(y(19)-y(39))...
    +(C(19,end)/A(end,end))*((-0.5/(u(end)-(u(end)^0.5)))-(A(end,1)*y(1)+A(end,2)*y(2)+A(end,3)*y(3)+A(end,4)*y(4)+A(end,5)*y(5)+A(end,6)*y(6)+A(end,7)*y(7)+A(end,8)*y(8)+A(end,9)*y(9)+A(end,10)*y(10)...
    +A(end,11)*y(11)+A(end,12)*y(12)+A(end,13)*y(13)+A(end,14)*y(14)+A(end,15)*y(15)+A(end,16)*y(16)+A(end,17)*y(17)+A(end,18)*y(18)+A(end,19)*y(19)+A(end,20)*y(20))));


dy(20)=(1/om)*(C(20,1)*y(1)+C(20,2)*y(2)+C(20,3)*y(3)+C(20,4)*y(4)+C(20,5)*y(5)+C(20,6)*y(6)+C(20,7)*y(7)+C(20,8)*y(8)+C(20,9)*y(9)+C(20,10)*y(10)...
  +C(20,11)*y(11)+C(20,12)*y(12)+C(20,13)*y(13)+C(20,14)*y(14)+C(20,15)*y(15)+C(20,16)*y(16)+C(20,17)*y(17)+C(20,18)*y(18)+C(20,19)*y(19)+C(20,20)*y(20)-la*(y(20)-y(40))...
    +(C(20,end)/A(end,end))*((-0.5/(u(end)-(u(end)^0.5)))-(A(end,1)*y(1)+A(end,2)*y(2)+A(end,3)*y(3)+A(end,4)*y(4)+A(end,5)*y(5)+A(end,6)*y(6)+A(end,7)*y(7)+A(end,8)*y(8)+A(end,9)*y(9)+A(end,10)*y(10)...
    +A(end,11)*y(11)+A(end,12)*y(12)+A(end,13)*y(13)+A(end,14)*y(14)+A(end,15)*y(15)+A(end,16)*y(16)+A(end,17)*y(17)+A(end,18)*y(18)+A(end,19)*y(19)+A(end,20)*y(20))));


dy(21)=(1/(1-om))*la*(y(1)-y(21));

dy(22)=(1/(1-om))*la*(y(2)-y(22));

dy(23)=(1/(1-om))*la*(y(3)-y(23));

dy(24)=(1/(1-om))*la*(y(4)-y(24));

dy(25)=(1/(1-om))*la*(y(5)-y(25));

dy(26)=(1/(1-om))*la*(y(6)-y(26));

dy(27)=(1/(1-om))*la*(y(7)-y(27));

dy(28)=(1/(1-om))*la*(y(8)-y(28));

dy(29)=(1/(1-om))*la*(y(9)-y(29));

dy(30)=(1/(1-om))*la*(y(10)-y(30));

dy(31)=(1/(1-om))*la*(y(11)-y(31));

dy(32)=(1/(1-om))*la*(y(12)-y(32));

dy(33)=(1/(1-om))*la*(y(13)-y(33));

dy(34)=(1/(1-om))*la*(y(14)-y(34));

dy(35)=(1/(1-om))*la*(y(15)-y(35));

dy(36)=(1/(1-om))*la*(y(16)-y(36));

dy(37)=(1/(1-om))*la*(y(17)-y(37));

dy(38)=(1/(1-om))*la*(y(18)-y(38));

dy(39)=(1/(1-om))*la*(y(19)-y(39));

dy(40)=(1/(1-om))*la*(y(20)-y(40));


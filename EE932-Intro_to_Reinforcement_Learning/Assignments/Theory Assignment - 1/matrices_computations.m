% computing the matrix crunching
D_a_1 = [1 3; 5 7]
b_a_1 = [17; 2]
I = eye(2)

D_a_2 = [7 13; 5 3]
b_a_2 = [2; 1]

x_7=[2;1]
theta_a_1 = inv(D_a_1'*D_a_1 + I)*D_a_1'*b_a_1
disp(theta_a_1)

theta_a_2 = inv(D_a_2'*D_a_2 + I)*D_a_2'*b_a_2
disp(theta_a_2)

mu_a_1 = theta_a_1'*x_7
mu_a_2 = theta_a_2'*x_7

disp(['mu_a_1:' + mu_a_1])
disp(['mu_a_2:' + mu_a_2])

% UCB computation
explore_a_1 = sqrt(x_7'*(inv(D_a_1'*D_a_1 + I))*x_7)
explore_a_2 = sqrt(x_7'*(inv(D_a_2'*D_a_2 + I))*x_7)
ucb_a_1= ( x_7'*theta_a_1 + sqrt(x_7'*(inv(D_a_1'*D_a_1 + I))*x_7) )
ucb_a_2= (x_7'*theta_a_2 + sqrt(x_7'*(inv(D_a_2'*D_a_2 + I))*x_7))
disp(['UCB of a_1:'+ ucb_a_1])
disp(['UCB of a_2:'+ ucb_a_2])

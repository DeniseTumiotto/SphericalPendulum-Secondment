q = rand(3,1);

q_skw = [0 q(3) -q(2); -q(3) 0 q(1); q(2) -q(1) 0];

disp(det(q_skw))

% x = q_skw\q;
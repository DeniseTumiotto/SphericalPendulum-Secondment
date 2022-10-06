[sols, params] = readUpto(10);

% evaluate error

error = evalErr(sols, params);
plotErr(error)

% animationSphere(sols(1), params(1));
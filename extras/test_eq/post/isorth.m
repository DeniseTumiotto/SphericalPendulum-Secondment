function rslt = isorth(R)

I_1 = R * R.';
I_2 = R.' * R;

tol = 1e-14;

if all(all(abs(I_1 - eye(3)) < tol)) && all(all(abs(I_2 - eye(3)) < tol))
    rslt = 1;
else
    rslt = 0;
end

end
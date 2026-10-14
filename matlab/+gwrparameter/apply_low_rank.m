function result = apply_low_rank(factors, values)
%APPLY_LOW_RANK Wendet A als P*(P.'*(W*v)) ohne dichte Matrix an.

values = values(:);
if numel(values) ~= factors.n+1
    error('gwrparameter:InvalidVectorSize', ...
        'values muss n+1 Eintraege besitzen.');
end

result = factors.factor*( ...
    factors.factor.'*(factors.weights.*values));

end

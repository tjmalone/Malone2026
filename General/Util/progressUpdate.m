function progressUpdate(numFiles)

    % Define count as a persistent variable
    persistent count;

    % Initialize count the first time the function is called
    if isempty(count)
        count = 0;
    end

    % Increment count
    count = count + 1;

    % print progress every 1%
    if mod(count, ceil(numFiles/100)) == 0
        fprintf('Progress: %.1f%%\n', (count / numFiles) * 100);
    end

end
function out = createSessionTable(measures, trialinfo, channels)
    out = [];
    % create metadata table based on Fieldtrip's trialinfo matrix
    meta = array2table(trialinfo, 'VariableNames', {'trial', 'condition', 'accuracy'});
    for measure = fieldnames(measures)'
        % average out the multiple estimates now that the measures have
        % been calculated. This is a no-op when not using multitapers.
        values = reshape(measures.(measure{1}), [], size(trialinfo, 1), size(channels, 2));
        values = squeeze(mean(values, 1));
        % For each calculated measure, create a table with columns for each
        % electrode...
        wide = array2table(values, 'VariableNames', channels);
        % ... & an extra column naming the current measure
        wide.measure = repelem(measure, size(wide, 1))';
        % then merge with the metadata table and convert the result from
        % wide to long form.
        long = stack([meta wide], channels, 'NewDataVariableName', 'value', 'IndexVariableName', 'channel');
        % ... and concatenate each measure-specific table to the output
        out = [out; long];
    end
end

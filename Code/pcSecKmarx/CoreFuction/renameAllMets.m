function model = renameAllMets(model)
% Normalize metabolite IDs and names by enforcing compartment

    % Construct mapping: compartment abbreviation → full name
    abbr2full = containers.Map(model.comps, model.compNames);

    comps = model.comps(:);
    compNames = model.compNames(:);

    for i = 1:numel(model.mets)
        met_id   = strtrim(model.mets{i});
        met_name = strtrim(model.metNames{i});
        met_id   = regexprep(met_id, '''+$', '');   
        met_name = regexprep(met_name, '''+$', ''); 


        % Extract compartment abbreviation from met ID
        tok_id = regexp(met_id, '\[([^\[\]]*)\]$', 'tokens', 'once');
        if isempty(tok_id)
            abbr = 'c';
            model.mets{i} = [met_id '[c]'];
        elseif isempty(tok_id{1}) || all(isspace(tok_id{1})) || all(tok_id{1} == 0)
            abbr = 'c';
            model.mets{i} = regexprep(met_id, '\[([^\[\]]*)\]$', '[c]');
        else
            abbr = tok_id{1};
        end

        % Skip entries where bracketed part is not a valid compartment
        if ~ismember(abbr, comps)
            continue;
        end
        fullName = abbr2full(abbr);

        % Process metNames compartment annotation
        tok_name_tail = regexp(met_name, '\[([^\[\]]+)\]$', 'tokens', 'once');

        if ~isempty(tok_name_tail)
            tail = strtrim(tok_name_tail{1});
            if isempty(tail)
                baseName = regexprep(met_name, '\[([^\[\]]+)\]$', '');
                model.metNames{i} = [baseName '[' fullName ']'];
                continue;
            end

            if strcmp(tail, fullName)
                continue;
            end

            % Remove existing valid compartment tag (abbr or full) and replace
            if ismember(tail, comps) || ismember(tail, compNames)
                baseName = regexprep(met_name, '\[([^\[\]]+)\]$', '');
                model.metNames{i} = [baseName '[' fullName ']'];
                continue;
            end

            % Tail exists but is not a valid compartment （append full name）
            model.metNames{i} = [met_name '[' fullName ']'];
        else
            % No compartment suffix present （append full name）
            model.metNames{i} = [met_name '[' fullName ']'];
        end
    end
end

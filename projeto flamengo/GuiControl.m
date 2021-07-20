classdef GuiControl < handle
    properties (Access = private)
        m_master_list = []; % list with users and passwords
        m_profile_quests = []; % profile questions
        m_app = []; % app
        m_user_name_flag = false; % flag for user identification
        m_user_password_flag = false; % flag for user identification
        m_autorization = false; %flag for valid authorization
        m_user_name = '' % user name;
        m_dialog = {};
        m_draw = [];
        m_password = ''; 
        m_server_name = 'Baica';
    end
    methods (Access = public)
        function this = GuiControl(app)
            this.m_app = app;
            this.m_master_list = readtable('Userid_MasterList.txt','Format','%s%s');
            this.Init();
            this.m_draw = Draw(app.ax_graph, app.table);
        end
        function res = IsPasswordAllowed(this)
            res = this.m_user_name_flag;
        end
        function CheckValidUserInput(this, user_input)
            if (~this.m_user_name_flag)
                this.m_user_name_flag = true;
                this.m_user_name = user_input;
                this.m_dialog(end-1) = strcat(this.m_dialog(end-1),{' '},user_input); % minus one due to space
                this.UpdateChat('Entre a sua senha:');
                return;
            elseif (~this.m_user_password_flag)
                this.m_dialog(end-1) = strcat(this.m_dialog(end-1),{' '},user_input);
                all_names = this.m_master_list{:,1};
                all_passwords = this.m_master_list{:,2};
                for i = 1:length(all_names)
                    name = all_names{i};
                    if (strcmp(name,this.m_user_name))
                        for j = 1:size(all_passwords)
                            password = all_passwords{j};
                            if (strcmp(password , this.m_password))
                                %User autorized
                                this.m_user_password_flag = true;
                                this.m_autorization = true;
                                return;
                            end
                        end
                    end
                end
            end
            % Try again
            this.UpdateChat('Autorização negada. Tente Novamente');
            this.Init();
        end
        function res = IsUserAutorized(this)
            res = this.m_autorization;
        end
        function Init(this)
            this.m_user_name_flag = false;
            this.m_user_name = '';
            this.m_password = '';
            this.UpdateChat('Entre seu userid:');
        end
        function Reset(this)
            this.m_user_name_flag = false;
            this.m_user_password_flag = false;
            this.m_dialog = {};
        end
        function InitDialogMode(this)
            this.Reset();
            if hour(datetime)<12
                Intro_Welcome='Bom dia,';
                
            elseif hour(datetime)<18
                Intro_Welcome='Boa Tarde,';
            else
                Intro_Welcome='Boa Noite,';
            end            
            Welcome = strcat(Intro_Welcome,this.m_user_name,'.',' Bem vindo ao Baica! Estou carregando o seu perfil e os últimos dados disponíveis. Qual é a sua primeira pergunta?');
            this.UpdateChat(Welcome);
            this.m_app.PerguntaDropDown.Value = {' '};
            this.m_profile_quests = readtable(strcat('Profile',this.m_user_name,'.txt'));
            this.m_app.PerguntaDropDown.Items = unique(this.m_profile_quests.Perguntas);
        end
        function UpdateChat(this,message,speaker)
            if nargin == 2
                speaker = this.m_server_name;
            end
            date_and_time = datetime;
            date_and_time = datetime(date_and_time,'Format','dd/MM/yyyy HH:mm');
            this.m_dialog(end+1) = strcat(char(date_and_time),{' '},speaker,{':'},{' '},message);
            this.m_dialog(end+1) = {' '};
            this.m_app.CHATcomoBAICATextArea.Value =  this.m_dialog;
            drawnow;
            scroll(this.m_app.CHATcomoBAICATextArea, 'bottom');
        end
        function ShowAnswer(this, question)
            res = 'error';
            TO = 1;
            while strcmp(res,'error')
                try
                    URL = 'https://v97vglvwj9.execute-api.us-east-1.amazonaws.com/dev/nlp/';
                    data= '{"text": "'+question+'" ,"user_id": "nome_identificador_do_usuario","debug_mode":"true","graphs":"true","app":"Baica"}';
                    data = char(data);
                    options=weboptions('Timeout', TO, 'RequestMethod', 'post', 'MediaType', 'application/json', 'HeaderFields', {'x-api-key' 'uKFINLIaQC8RMGGAtDk4j7xpB7o5pnaK3i38fPZq'
                        });
                    options.ArrayFormat='json';
                    res = webwrite(URL,data,options);
                catch ME
                    ME.message;
                    ME.stack;
                    ME.identifier;
                    if strcmp(ME.identifier,'MATLAB:webservices:Timeout')
                        res = 'error';
                        TO = 1;
                    end
                end
            end
            this.UpdateChat(res.answer);
            this.m_draw.Reset();
            % check if exists messasges to draw
            fn = fieldnames(res);
            if (ismember('graphs',fn))
                this.m_draw.draw(res.graphs);
            end
        end
        function BuildPassowrd(this,letter)
            this.m_password(end+1) = letter;
        end
        function res = GetUserName(this)
            res = this.m_user_name;
        end
    end
    
end
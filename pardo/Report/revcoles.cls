%%
%% This is file `revcoles.cls', REVCOLES DOCUMENT CLASS --
%%
%% Designed by Prof. CAMILO CUBIDES GARZÓN, eccubidesg@unal.edu.co, eccubidesg@gmail.com
%%
%% Copyright UNIVERSIDAD NACIONAL DE COLOMBIA® - NATIONAL UNIVERSITY OF COLOMBIA®
%%
%%           REVISTA COLOMBIANA DE ESTADÍSTICA - COLOMBIAN STATISTICS JOURNAL
%%
%% This file was generated from file(s) of the LaTeX base system.
%% --------------------------------------------------------------
%%
%% \CharacterTable
%%  {Upper-case    \A\B\C\D\E\F\G\H\I\J\K\L\M\N\O\P\Q\R\S\T\U\V\W\X\Y\Z
%%   Lower-case    \a\b\c\d\e\f\g\h\i\j\k\l\m\n\o\p\q\r\s\t\u\v\w\x\y\z
%%   Digits        \0\1\2\3\4\5\6\7\8\9
%%   Exclamation   \!     Double quote  \"     Hash (number) \#
%%   Dollar        \$     Percent       \%     Ampersand     \&
%%   Acute accent  \'     Left paren    \(     Right paren   \)
%%   Asterisk      \*     Plus          \+     Comma         \,
%%   Minus         \-     Point         \.     Solidus       \/edition
%%   Colon         \:     Semicolon     \;     Less than     \<
%%   Equals        \=     Greater than  \>     Question mark \?
%%   Commercial at \@     Left bracket  \[     Backslash     edition\\
%%   Right bracket \]     Circumflex    \^     Underscore    \_
%%   Grave accent  \`     Left brace    \{     Vertical bar  \|
%%   Right brace   \}     Tilde         \~}
\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{revcoles}[2009/04/22 v2.5 document class for multi authored article for REVCOLES journal, Symposium of statistic, projects of pregrade in statistic, and report in a course at the National University of Colombia]

\RequirePackage[T1]{fontenc}                        % fuente utilizada para incorporar simbolos acentuados y otros símbolos
\RequirePackage{amsmath,amsthm,amssymb,latexsym}    % herramientas para manejo matemático y símbolos matemáticos especiales
             
\RequirePackage{verbatim}                           % herramienta para incluir texto en fuente mono-espaciada y para comentarios
\RequirePackage{enumerate}                          % herramienta para el manejo de listas enumeradas y encajadas
\RequirePackage{array}                              % utilizado para darle formato al encabezado de las tablas
\RequirePackage{multicol}
\RequirePackage{ifpdf}
\ifpdf
    \RequirePackage{graphicx,color}%                % herramientas para inclusión de gráficos y el manejo del color
    \DeclareGraphicsExtensions{.pdf,.png,.jpg}
\else
    \RequirePackage[dvips]{graphicx,color}%
    \DeclareGraphicsExtensions{.eps}
\fi

\def\BibTeX{B{\footnotesize IB}\TeX}            % para agregar la palabra BibTeX en estilo LaTeX
\def\MiKTeX{M{\footnotesize iK}\TeX}            % para agregar la palabra MikTeX en estilo LaTeX
\def\pdfLaTeX{\textsf{pdf}\LaTeX}               % para agregar la palabra pdfLaTeX en estilo LaTeX
\def\emphred#1{{\color{red}#1}}                 % para resaltar texto dudoso a corregir

\definecolor{lightgray}{gray}{0.5}

\def\Hy@colorlink{}
\def\Hy@endcolorlink{}

\def\@labeljournal{}
\def\@firstpage{1}
\def\@lastpage{0}
\def\@volume{0}
\def\@number{0}
\def\@month{1}
\def\@monthaux{0}
\def\@year{0000}
\def\@monthreceived{0}
\def\@monthaccepted{0}
\def\@yearreceived{0000}
\def\@yearaccepted{0000}
\def\@subdivision{}
\def\@division{}
\def\@institution{}
\def\@city{}
\def\@country{}
\def\@code{}
\def\@maintitle{}
\def\@secondtitle{}
\def\@shorttitle{}
\def\@symposiumname{}
\def\@symposiumnumber{}
\def\@symposiumcity{}
\def\@symposiumcountry{}
\def\@symposiumdate{}
\def\@labelproject{}
\def\@labelUniversity{}
\def\@labelFaculty{}
\def\@labelDepartment{}
\def\@labelProgram{}
\def\@labelSubject{}
\def\@labelCity{Bogotá}
\def\@labelProject{}
\newcounter{@symposium}
\newcounter{@instnumber}
\newcounter{@inst}
\newcounter{@last}
\def\@inst{}
\newcounter{@author}
\setcounter{@author}{1}
\def\@firstname{}
\def\@surname{}
\def\@numberinstitution{}
\def\@affiliation{}
\def\@email{}
\def\@names{}
\def\@namesaux{}
\def\@namestablecontents{}
\def\@namestablecontentsaux{}
\xdef\labeleditorial{}
\def\footmark{}
\def\keywordsname{}
\def\principalabstractname{}
\def\secondaryabstractname{}
\def\principalkeywordsname{}
\def\secondarykeywordsname{}
\def\¿{?`}
\newif\if@report\@reportfalse
\newif\if@project\@projectfalse
\newif\if@referees\@refereesfalse
\newif\if@receivedaccepted\@receivedacceptedfalse
\newif\if@shortauthor\@shortauthorfalse
\newif\if@languagespanish\@languagespanishfalse
\newif\if@languageenglish\@languageenglishfalse
\newif\if@languageportuguese\@languageportuguesefalse

\def\@prop{Proposición}
\def\@teor{Teorema}
\def\@coro{Corolario}
\def\@resu{Resultado}
\def\@lema{Lema}
\def\@defi{Definición}
\def\@exam{Ejemplo}
\def\@exams{Ejemplos}
\def\@exer{Ejercicio}
\def\@exers{Ejercicios}
\def\@note{Nota}
\def\@nota{Notación}
\def\@obse{Observación}
\def\@cond{Condición}

\theoremstyle{plain}
    \newtheorem{prop}{\@prop}%
    \newtheorem{teor}{\@teor}%
    \newtheorem{coro}{\@coro}%
    \newtheorem{lema}{\@lema}%
    \newtheorem{resu}{\@resu}%
\theoremstyle{definition}
    \newtheorem{defi}{\@defi}%
    \newtheorem{exam}{\@exam}%
    \newtheorem{exams}{\@exams}%
    \newtheorem{exer}{\@exer}%
    \newtheorem{exers}{\@exers}%
\theoremstyle{remark}
    \newtheorem{note}{\@note}%
    \newtheorem{nota}{\@nota}%
    \newtheorem{obse}{\@obse}%
    \newtheorem{cond}{\@cond}%

\newcommand\@ptsize{}
\newif\if@restonecol
\newif\if@titlepage
\@titlepagefalse

\def\@smallpage{
    \input{size1\@ptsize.clo}
    \setlength\paperheight{11in}%{240mm}
    \setlength\paperwidth{8.5in}%{165mm}
    \setlength\textheight{197mm}
    \setlength\textwidth{125mm}
    \setlength\parskip{1mm}
    \setlength\marginparwidth{2mm}
    \setlength\headsep{6mm}
    \setlength\oddsidemargin{23mm}
    \setlength\evensidemargin{23mm}
    \setlength\voffset{-1mm}
}

\def\@bigpage{
    \setlength\paperheight{11in}%
    \setlength\paperwidth{8.5in}%
    \setlength\textheight{225mm}
    \setlength\textwidth{160mm}
    \setlength\topmargin{0cm}
    \setlength\parskip{2mm}
    \setlength\marginparwidth{2mm}
    \setlength\headsep{6mm}
    \setlength\oddsidemargin{11mm}
    \setlength\evensidemargin{-4.5mm}
    \setlength\voffset{-5.4mm}
}

\if@compatibility
  \renewcommand\@ptsize{0}
\else
\DeclareOption{10pt}{\renewcommand\@ptsize{0}}
\fi
\DeclareOption{11pt}{\renewcommand\@ptsize{1}}
\DeclareOption{12pt}{\renewcommand\@ptsize{2}}

\if@compatibility\else
\DeclareOption{revcolespaper}{
    \@smallpage
}
\fi

\if@compatibility\else
\DeclareOption{oneside}{\@twosidefalse \@mparswitchfalse}
\fi
\DeclareOption{twoside}{\@twosidetrue  \@mparswitchtrue}
\DeclareOption{draft}{\setlength\overfullrule{5pt}}

\if@compatibility\else
\DeclareOption{final}{\setlength\overfullrule{0pt}}
\fi
\DeclareOption{notitlepage}{\@titlepagefalse}
\DeclareOption{twocolumn}{\@twocolumntrue}
\DeclareOption{leqno}{\input{leqno.clo}}
\DeclareOption{fleqn}{\input{fleqn.clo}}
\DeclareOption{openbib}{%
  \AtEndOfPackage{%
   \renewcommand\@openbib@code{%
      \advance\leftmargin\bibindent
      \itemindent -\bibindent
      \listparindent \itemindent
      \parsep \z@
      }%
   \renewcommand\newblock{\par}}%
}
%\if@compatibility\else
\DeclareOption{spanish}{
    \@languagespanishtrue
    \@languageenglishfalse
    \@languageportuguesefalse
    \def\principalabstractname{Resumen}
    \def\secondaryabstractname{Abstract}
    \def\principalkeywordsname{Palabras clave}
    \def\secondarykeywordsname{Key words}
    \def\@prop{Proposición}
    \def\@teor{Teorema}
    \def\@coro{Corolario}
    \def\@resu{Resultado}
    \def\@lema{Lema}
    \def\@defi{Definición}
    \def\@defis{Definiciones}
    \def\@exam{Ejemplo}
    \def\@exams{Ejemplos}
    \def\@exer{Ejercicio}
    \def\@exers{Ejercicios}
    \def\@note{\bf{\emph{Nota}}}
    \def\@nota{\bf{\emph{Notación}}}
    \def\@obse{\bf{\emph{Observación}}}
    \def\@cond{\bf{\emph{Condición}}}
    \def\@appendix{Apéndice}
}
\DeclareOption{english}{
    \@languagespanishfalse
    \@languageenglishtrue
    \@languageportuguesefalse
    \def\principalabstractname{Abstract}
    \def\secondaryabstractname{Resumen}
    \def\principalkeywordsname{Key words}
    \def\secondarykeywordsname{Palabras clave}
    \def\@prop{Proposition}
    \def\@teor{Theorem}
    \def\@coro{Corollary}
    \def\@lema{Lema}
    \def\@defi{Definition}
    \def\@exam{Example}
    \def\@exer{Exercise}
    \def\@appendix{Appendix}
}
\DeclareOption{portuguese}{
    \@languagespanishfalse
    \@languageenglishfalse
    \@languageportuguesetrue
    \def\principalabstractname{Resumo}
    \def\secondaryabstractname{Abstract}
    \def\principalkeywordsname{Palavras chave}
    \def\secondarykeywordsname{Key words}
    \def\@prop{Proposição}
    \def\@teor{Teorema}
    \def\@coro{Corolário}
    \def\@lema{Lema}
    \def\@defi{Definição}
    \def\@exam{Exemplo}
    \def\@exer{Exercício}
    \def\@appendix{Apêndice}
}
%\fi

\if@compatibility\else
\DeclareOption{author}{\def\@labeljournal{\labelauthor}}
\DeclareOption{referees}{\def\@labeljournal{\labelreferees}\@refereestrue}
\DeclareOption{edition}{\def\@labeljournal{\labeledition}\@receivedacceptedtrue}
\DeclareOption{editorial}{\xdef\@labeljournal{\labeleditorial}}
\DeclareOption{publication}{\def\@labeljournal{\labelpublication}\@receivedacceptedtrue}
\DeclareOption{symposium}{\@bigpage\def\@labeljournal{\labelsymposium}}
\DeclareOption{project}{\@bigpage\def\@labeljournal{\labelproject}\@projecttrue}
\DeclareOption{report}{\@bigpage\def\@labeljournal{\labelreport}\@reporttrue}
\fi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\ExecuteOptions{10pt,revcolespaper,twoside,onecolumn,final,author,spanish}%
\ProcessOptions


\def\layoutmultilanguageitems{
    \def\labelitemi{\textbullet}%
    \def\labelitemii{$\m@th\diamond$}%
    \def\labelitemiii{\leavevmode\hbox to 1.2ex {\hss\vrule height .9ex width .7ex depth -.2ex\hss}}%
    \def\labelitemiv{$\m@th\circ$}%
    \def\labelenumi{\arabic{enumi}.}%
    \def\labelenumii{\alph{enumii})}%
    \def\labelenumiii{\roman{enumiii}.}%
    \def\labelenumiv{\Alph{enumiv})}%
}

\if@languageenglish
    \usepackage[english]{babel}
    \def\@seccntformat#1{\csname the#1\endcsname\space}
    \addto\captionsenglish{
        \def\proofname{\bf{\emph{Proof}}}
    }
    \layoutmultilanguageitems
    \let\@afterindentfalse\@afterindenttrue
\else\if@languagespanish
        \usepackage[spanish]{babel}
        \addto\layoutspanish{\def\@seccntformat#1{\csname the#1\endcsname.\space}\layoutmultilanguageitems}
        \decimalpoint
        \deactivatetilden
        \addto\captionsspanish{
            \def\proofname{\bf{\emph{Demostración}}}
            \def\tablename{Tabla}
        }
\else\if@languageportuguese
        \usepackage[portuguese]{babel}
        \def\@seccntformat#1{\csname the#1\endcsname.\space}
        \addto\captionsportuguese{
            \def\proofname{\bf{\emph{Demonstração}}}
        }
        \layoutmultilanguageitems
        \let\@afterindentfalse\@afterindenttrue
\fi
\fi
\fi



\def\appendix{\par
    \setcounter{section}{0}%
    \setcounter{subsection}{0}%
    \if@languagespanish
    	\def\thesection {\@appendix\ \@Alph\c@section}
    	\def\thesubsection   {\@Alph\c@section.\@arabic\c@subsection}
	\else
    	\def\thesection {\@appendix\ \@Alph\c@section.}
   		\def\thesubsection   {\@Alph\c@section.\@arabic\c@subsection.}
	\fi
}



\AtBeginDocument{\ifx\@labeljournal\labeleditorial\pagestyle{empty}\else\pagestyle{revcoleshead}\thispagestyle{plain}\fi\def\@roman#1{\romannumeral #1}}
\AtEndDocument{\label{@lastpage}\setcounter{@last}{\value{page}}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\def\labeldocument{\@ifnextchar[\@valueDocument{}}
\def\@valueDocument[#1]{\setkeys{values}{#1}}

\define@key{values}{firstpage}{\def\@firstpage{#1}\setcounter{page}{\@firstpage}}
\define@key{values}{volume}{\def\@volume{#1}}
\define@key{values}{number}{\def\@number{#1}}
\define@key{values}{month}{\def\@month{#1}\def\@monthaux{#1}}
\define@key{values}{year}{\def\@year{#1}}
\define@key{values}{monthreceived}{\def\@monthaux{#1}\xdef\@monthreceived{\month}}
\define@key{values}{monthaccepted}{\def\@monthaux{#1}\xdef\@monthaccepted{\month}}
\define@key{values}{yearreceived}{\def\@yearreceived{#1}}
\define@key{values}{yearaccepted}{\def\@yearaccepted{#1}}
\define@key{values}{symposiumtitle}{\def\@symposiumtitle{#1}}
\define@key{values}{symposiumname}{\if#1\empty\else\def\@symposiumname{\\#1}\fi}
\define@key{values}{symposiumnumber}{\setcounter{@symposium}{#1}\def\@symposiumnumber{\Roman{@symposium}}}%
\define@key{values}{symposiumcity}{\def\@symposiumcity{#1}}
\define@key{values}{symposiumcountry}{\def\@symposiumcountry{#1}}
\define@key{values}{symposiumdate}{\def\@symposiumdate{#1}}
\define@key{values}{label}{\def\@labelProject{\\#1}}
\define@key{values}{university}{\def\@labelUniversity{#1\\}}
\define@key{values}{faculty}{\def\@labelFaculty{#1\\}}
\define@key{values}{department}{\def\@labelDepartment{#1\\}}
\define@key{values}{program}{\def\@labelProgram{#1\\}}
\define@key{values}{subject}{\def\@labelSubject{#1}}
\define@key{values}{city}{\def\@labelCity{#1}}


\def\Monthpublication{\ifcase\@month\or Enero\or Febrero\or Marzo\or Abril\or Mayo\or Junio\or Julio\or Agosto\or Septiembre\or Octubre\or Noviembre\or Diciembre\fi}
\def\month{\ifcase\@monthaux\or enero\or febrero\or marzo\or abril\or mayo\or junio\or julio\or agosto\or septiembre\or octubre\or noviembre\or diciembre\fi}


\def\title{\@ifnextchar[\@valueTitle{}}%
\def\@valueTitle[#1]{%
    \setkeys{valuesTitle}{#1}
    \@labeljournal
    \begin{center}\Large\bfseries \@maintitle \end{center}
    \begin{center}\normalsize\bfseries \@secondtitle \end{center}
    \ifx\@empty\@shorttitle
        \gdef\rightmark{\@maintitle}
    \else
        \gdef\rightmark{\@shorttitle}
    \fi
}

\def\infojournal{\footnotesize%
\hspace*{-4cm}\mbox{\string\textbf\{\@namestablecontents\}\textbackslash\textbackslash}\newline%
\hspace*{-4cm}\mbox{\string\emph\{\@maintitle\}\%}\newline%
\hspace*{-4cm}\mbox{\%\string\dotfill\textbackslash\textbackslash\string\phantom\{\string\hspace\{-0.5ex\}\}\%}\newline%
\hspace*{-4cm}\mbox{\string\dotfill\space\@firstpage-\pageref{@lastpage}\textbackslash\textbackslash[-1mm]}
}

\define@key{valuesTitle}{maintitle}{\def\@maintitle{#1}}
\define@key{valuesTitle}{secondtitle}{\def\@secondtitle{#1}}
\define@key{valuesTitle}{shorttitle}{\def\@shorttitle{#1}}

\newenvironment{authors}[1][]{%
    \ifx#1\@empty
    \else
        \xdef\leftmark{#1}
        \@shortauthortrue
    \fi
    \if@project\begin{multicols}{2}\fi
    \if@referees\else\setcounter{footnote}{0}\def\thefootnote{\@alph\c@footnote}\begin{center}\scshape\fi
}{
    \if@referees
    \else
        \end{center}\setcounter{footnote}{0}\def\thefootnote{\@arabic\c@footnote}
        \if@shortauthor\else\xdef\leftmark{\@names}\fi
    \fi
    \if@project\end{multicols}\fi
    %\infojournal%para crear el índice de autores contraportada
}

\def\author{\@ifnextchar[\@valueAuthor{}}
\def\@valueAuthor[#1]{%
    \setkeys{valuesAuthor}{#1}\@name%
    \def\@firstname{}%
    \def\@surname{}%
    \def\@numberinstitution{}%
    \def\@affiliation{}%
    \def\@email{}%
	\def\@code{}%
}

\define@key{valuesAuthor}{firstname}{\def\@firstname{#1}}
\define@key{valuesAuthor}{surname}{\def\@surname{#1}}
\define@key{valuesAuthor}{numberinstitution}{\def\@numberinstitution{#1,}\stepcounter{@instnumber}}
\define@key{valuesAuthor}{affiliation}{\def\@affiliation{#1}}
\define@key{valuesAuthor}{email}{\def\@email{.\ E-mail: #1}}
\define@key{valuesAuthor}{code}{\def\@code{#1}}

\def\@name{
\if@report
    \ifnum\value{@author} = 1%
        \mbox{\@firstname\ \@surname}\nolinebreak\unskip\\%\footnote{\@affiliation\@code\@email}%
            \xdef\@names{\@firstname\ \@surname}%
            \xdef\@namesaux{\@firstname\ \@surname}%
        \else%
            \unskip\unskip, \mbox{\@firstname\ \@surname}\nolinebreak\unskip\unskip\footnote{\@affiliation\@code\@email}%
            \xdef\@names{\@namesaux\ \& \@firstname\ \@surname}%
            \xdef\@namesaux{\@namesaux, \@firstname\ \@surname}%
    \fi%
    \stepcounter{@author}%
\else
\if@project
    \ifnum\value{@author} = 1%
        \mbox{\@firstname\ \@surname}\nolinebreak\unskip\\%\footnote{\@affiliation \@email}\\%
        {\small \if@languagespanish Autor \else Author \fi}\\
        \xdef\@names{\@firstname\ \@surname}%
        \xdef\@namesaux{\@firstname\ \@surname}%
    \else%
        \mbox{\@firstname\ \@surname}\nolinebreak\unskip\\%\footnote{\@affiliation \@email}\\%
        {\small \if@languagespanish Director \else Advisor \fi}
    \fi%
    \stepcounter{@author}%
\else
\if@referees
    \gdef\leftmark{\textit{\rightmark}}%
\else
    \ifnum\value{@author} = 1%
        \mbox{\@firstname\space\@surname}\nolinebreak\unskip\\%$^{\@numberinstitution}$\nolinebreak\unskip\footnote{\@affiliation \@email}%
        \xdef\@names{\@firstname\ \@surname}%
        \xdef\@namesaux{\@firstname\ \@surname}%
        \xdef\@namestablecontents{\@firstname\ \@surname}%
        \xdef\@namestablecontentsaux{\@firstname\ \@surname}%
    \else%
        \unskip\unskip, \mbox{\@firstname\space\@surname}\\%\nolinebreak\unskip$^{\@numberinstitution}$\nolinebreak\unskip\footnote{\@affiliation \@email}%
        \xdef\@names{\@namesaux\ \& \@firstname\ \@surname}%
        \xdef\@namesaux{\@namesaux, \@firstname\ \@surname}%
        \xdef\@namestablecontents{\@namestablecontentsaux\ \textbackslash\& \@firstname\ \@surname}%
        \xdef\@namestablecontentsaux{\@namestablecontentsaux, \@firstname\ \@surname}%
    \fi%
    \stepcounter{@author}%
\fi%
\fi%
\fi%
}

\newenvironment{institutions}{%
\if@referees\else\relax\setcounter{@inst}{0}\fi
\begin{center}\footnotesize\scshape
}{
\end{center}\separator
}

\def\institute{\@ifnextchar[\@valueInstitute{}}
\def\@valueInstitute[#1]{
\setkeys{valuesInst}{#1}%
\if@referees
\else
    \ifnum \value{@instnumber} = 0
        \@subdivision \@division \@institution \@city \@country\\
    \else
        \stepcounter{@inst}%\enspace%\unskip%\ignorespaces
        \noindent$^\the@inst${\ignorespaces\@subdivision \@division \@institution \@city \@country}\\[1ex]
    \fi
\fi%
\def\@subdivision{}%
\def\@division{}%
\def\@institution{}%
\def\@city{}%
\def\@country{}%
}

\define@key{valuesInst}{subdivision}{\def\@subdivision{#1, }}
\define@key{valuesInst}{division}{\def\@division{#1, }}
\define@key{valuesInst}{institution}{\def\@institution{#1, }}
\define@key{valuesInst}{city}{\def\@city{#1, }}
\define@key{valuesInst}{country}{\def\@country{#1}}

\newcommand\separator{\vspace{-4pt}\centerline{\color{lightgray}\rule{45mm}{0.5pt}}\vspace{13pt}}

\def\multititle#1{\multicolumn{1}{c}{\renewcommand{\tabcolsep}{0pt}\begin{tabular}{c}#1\end{tabular}}}
\def\multititlerot#1{\multicolumn{1}{c}{\renewcommand{\tabcolsep}{0pt}\rotatebox{90}{\begin{tabular}{l}#1\end{tabular}}}}
\def\multititlefinal#1{\multicolumn{1}{c}{\renewcommand{\tabcolsep}{0pt}\begin{tabular}{c}#1\end{tabular}}}
\def\multititles#1#2{\multicolumn{#1}{c}{\renewcommand{\tabcolsep}{0pt}\begin{tabular}{c}#2\end{tabular}}}
\def\multititlefinals#1#2{\multicolumn{#1}{c}{\renewcommand{\tabcolsep}{0pt}\begin{tabular}{c}#2\end{tabular}}}
\def\multititlel#1{\multicolumn{1}{l}{\renewcommand{\tabcolsep}{0pt}\begin{tabular}{l}#1\end{tabular}}}
\def\multititlefinall#1{\multicolumn{1}{l}{\renewcommand{\tabcolsep}{0pt}\begin{tabular}{l}#1\end{tabular}}}

\newenvironment{mainabstract}{%
      \def\abstractname{\principalabstractname}
      \def\keywordsname{\principalkeywordsname}
        \small
        \begin{center}%
          {\bfseries \abstractname\vspace{-.5em}\vspace{\z@}}%
        \end{center}%
        \quotation
}{\endquotation}

\newenvironment{secondaryabstract}{%
      \def\abstractname{\secondaryabstractname}
      \def\keywordsname{\secondarykeywordsname}
        \small
        \begin{center}%
          {\bfseries \abstractname\vspace{-.5em}\vspace{\z@}}%
        \end{center}%
        \quotation
}{\endquotation}

\def\keywords#1{\\[1.5ex]{\bfseries\textit\keywordsname:\;\unskip}#1.}

\def\ps@revcoleshead{
      \def\@oddfoot{\hfil\footnotesize \footmark}%
      \def\@evenfoot{\hfil\footnotesize \footmark}%
      \def\@evenhead{\thepage\hfil{\small\leftmark}}%
      \def\@oddhead{{\small\textit\rightmark}\hfil\thepage}%
      \let\@mkboth\markboth
    \def\sectionmark##1{%
      \markboth {\MakeUppercase{%
        \ifnum \c@secnumdepth >\z@
          \thesection\quad
        \fi
        ##1}}{}}%
    \def\subsectionmark##1{%
      \markright {%
        \ifnum \c@secnumdepth >\@ne
          \thesubsection\quad
        \fi
        ##1}}
}

\newcommand{\receivedaccepted}[2]{%
    %\ \medskip\\%
    \begin{center}
        \textsf{\big[Recibido: #1 --- Aceptado: #2\big]}
    \end{center}
}

\def\references#1{\bibliography{#1}}

\def\bibliography#1{%
    \if@receivedaccepted
        \receivedaccepted{\@monthreceived\ de\ \@yearreceived}{\@monthaccepted\ de\ \@yearaccepted}
    \fi
    \bibliographystyle{agsm}%
    \if@filesw
        \immediate\write\@auxout{\string\bibdata{#1}}%
    \fi
    \@input@{\jobname.bbl}
}

\def\labelauthor{
\vspace*{-14mm}%
\begin{flushright}
    \small%
    \textbf{Revista Colombiana de Estadística}\\[1pt]
    Edición para autor(es)\\[-9pt]
    {\color{lightgray}\underline{\mbox{\hphantom{Edición para autor(es)}}}}
\end{flushright}\vspace{8mm}
    \def\footmark{Revista Colombiana de Estadística, Edición para autor(es)}
}

\def\labelreferees{
\vspace*{-14mm}%
\begin{flushright}
    \small%
    \textbf{Revista Colombiana de Estadística}\\[1pt]
    Edición para arbitraje\\[-9pt]
    {\color{lightgray}\underline{\mbox{\hphantom{Edición para arbitraje}}}}
\end{flushright}\vspace{8mm}
    \def\footmark{Revista Colombiana de Estadística, edición para arbitraje}
}

\def\labeledition{
\vspace*{-14mm}%
\begin{flushright}
    \small%
    \textbf{Revista Colombiana de Estadística}\\[1pt]
    Artículo en edición, pp.\ \@firstpage\ a \pageref{@lastpage}.\ \@year\\[-9pt]
    {\color{lightgray}\underline{\mbox{\hphantom{Artículo en edición, pp.\ \@firstpage\ a \pageref{@lastpage}.\ \@year}}}}
\end{flushright}\vspace{8mm}
\def\footmark{Revista Colombiana de Estadística. En edición (\@year) \@firstpage--\pageref{@lastpage}}
}

\def\labelpublication{
\vspace*{-14mm}%
\begin{flushright}
    \small%
    \textbf{Revista Colombiana de Estadística}\\[1pt]
    \Monthpublication\ \@year, volumen \@volume,\ no.\ \@number,\ pp.\ \@firstpage\ a \pageref{@lastpage}\\[-9pt]
    {\color{lightgray}\underline{\mbox{\hphantom{\Monthpublication\ \@year, volumen \@volume,\ no.\ \@number,\ pp.\ \@firstpage\ a \pageref{@lastpage}}}}}
\end{flushright}\vspace{8mm}
    \def\footmark{Revista Colombiana de Estadística {\bfseries\@volume} (\@year) \@firstpage--\pageref{@lastpage}}
}

\def\labelsymposium{
\vspace*{-14mm}%
\begin{flushright}
    \small%
    \textbf{\@symposiumnumber\space\@symposiumtitle\@symposiumname}\\[1pt]
    \@symposiumcity, \@symposiumcountry, \@symposiumdate\ de\ \@year\\[-9pt]
    {\color{lightgray}\underline{\mbox{\hphantom{\@symposiumcity, \@symposiumcountry, \@symposiumdate\ de\ \@year}}}}
\end{flushright}\vspace{8mm}
    \def\footmark{\@symposiumnumber\ Simposio de Estadística (\@year)}
}

\def\labelproject{
\vspace*{-14mm}%
\begin{flushright}
    \small%
    \textbf{University of Warwick\\
        Centre for Complexity Science\@labelProject}\\[1pt]
   \@year\\[-9pt]
    {\color{lightgray}\underline{\mbox{\hphantom{\@year}}}}
\end{flushright}\vspace{8mm}
    %\def\footmark{Trabajo de grado, Carrera de Estadística (\@year)}
}

\def\labelreport{
\vspace*{-14mm}%
\begin{flushright}
    \small%
    \textbf{\@labelUniversity\@labelFaculty\@labelDepartment\@labelProgram\@labelSubject}\\[1pt]
    %\@labelCity,\ \month\ de \@year\\[-9pt]
   % {\color{lightgray}\underline{\mbox{\hphantom{Bogotá, \month\ de \@year}}}}
\end{flushright}\vspace{8mm}
   % \def\footmark{\@labelSubject\ (\@year)}
}

\def\appendixsection#1{
    \section{}
    \ifx#1\@empty
    \else
        \vspace{-8pt}
        \noindent\textbf{\begin{Large}#1\end{Large}}\bigskip
    \fi
}

%\input{size1\@ptsize.clo}
\setlength\lineskip{1\p@}
\setlength\normallineskip{1\p@}
\renewcommand\baselinestretch{}
\setlength\parskip{3\p@ \@plus \p@}
\@lowpenalty   51
\@medpenalty  151
\@highpenalty 301
\setcounter{topnumber}{2}
\renewcommand\topfraction{.7}
\setcounter{bottomnumber}{1}
\renewcommand\bottomfraction{.3}
\setcounter{totalnumber}{3}
\renewcommand\textfraction{.2}
\renewcommand\floatpagefraction{.5}
\setcounter{dbltopnumber}{2}
\renewcommand\dbltopfraction{.7}
\renewcommand\dblfloatpagefraction{.5}

\setcounter{secnumdepth}{3}
\newcounter {section}
\newcounter {subsection}[section]
\newcounter {subsubsection}[subsection]
\newcounter {paragraph}[subsubsection]
\newcounter {subparagraph}[paragraph]

\if@languagespanish
    \renewcommand\thesection {\@arabic\c@section}
    \renewcommand\thesubsection   {\thesection.\@arabic\c@subsection}
    \renewcommand\thesubsubsection{\thesubsection.\@arabic\c@subsubsection}
    \renewcommand\theparagraph    {\thesubsubsection.\@arabic\c@paragraph}
    \renewcommand\thesubparagraph {\theparagraph.\@arabic\c@subparagraph}
\else
    \renewcommand\thesection {\@arabic\c@section.}
    \renewcommand\thesubsection   {\thesection\@arabic\c@subsection.}
    \renewcommand\thesubsubsection{\thesubsection\@arabic\c@subsubsection.}
    \renewcommand\theparagraph    {\thesubsubsection\@arabic\c@paragraph.}
    \renewcommand\thesubparagraph {\theparagraph\@arabic\c@subparagraph.}
\fi
\newcommand\section{\@startsection {section}{1}{\z@}%
                                   {-3.5ex \@plus -1ex \@minus -.2ex}%
                                   {2.3ex \@plus.2ex}%
                                   {\normalfont\Large\bfseries}}
\newcommand\subsection{\@startsection{subsection}{2}{\z@}%
                                     {-3.25ex\@plus -1ex \@minus -.2ex}%
                                     {1.5ex \@plus .2ex}%
                                     {\normalfont\large\bfseries}}
\newcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}%
                                     {-3.25ex\@plus -1ex \@minus -.2ex}%
                                     {1.5ex \@plus .2ex}%
                                     {\normalfont\normalsize\bfseries}}
\newcommand\paragraph{\@startsection{paragraph}{4}{\z@}%
                                    {3.25ex \@plus1ex \@minus.2ex}%
                                    {-1em}%
                                    {\normalfont\normalsize\bfseries}}
\newcommand\subparagraph{\@startsection{subparagraph}{5}{\parindent}%
                                       {3.25ex \@plus1ex \@minus .2ex}%
                                       {-1em}%
                                      {\normalfont\normalsize\bfseries}}

\if@twocolumn
  \setlength\leftmargini  {2em}
\else
  \setlength\leftmargini  {2.5em}
\fi
\leftmargin  \leftmargini
\setlength\leftmarginii  {2.2em}
\setlength\leftmarginiii {1.87em}
\setlength\leftmarginiv  {1.7em}
\if@twocolumn
  \setlength\leftmarginv  {.5em}
  \setlength\leftmarginvi {.5em}
\else
  \setlength\leftmarginv  {1em}
  \setlength\leftmarginvi {1em}
\fi
\setlength  \labelsep  {.5em}
\setlength  \labelwidth{\leftmargini}
\addtolength\labelwidth{-\labelsep}
\@beginparpenalty -\@lowpenalty
\@endparpenalty   -\@lowpenalty
\@itempenalty     -\@lowpenalty
\renewcommand\theenumi{\@arabic\c@enumi}
\renewcommand\theenumii{\@alph\c@enumii}
\renewcommand\theenumiii{\@roman\c@enumiii}
\renewcommand\theenumiv{\@Alph\c@enumiv}
%\newcommand\labelenumi{\theenumi.}
%\newcommand\labelenumii{(\theenumii)}
%\newcommand\labelenumiii{\theenumiii.}
%\newcommand\labelenumiv{\theenumiv.}
\renewcommand\p@enumii{\theenumi}
\renewcommand\p@enumiii{\theenumi(\theenumii)}
\renewcommand\p@enumiv{\p@enumiii\theenumiii}
%\newcommand\labelitemi{\textbullet}
%\newcommand\labelitemii{\normalfont\bfseries \textendash}
%\newcommand\labelitemiii{\textasteriskcentered}
%\newcommand\labelitemiv{\textperiodcentered}
\newenvironment{description}
               {\list{}{\labelwidth\z@ \itemindent-\leftmargin
                        \let\makelabel\descriptionlabel}}
               {\endlist}
\newcommand*\descriptionlabel[1]{\hspace\labelsep
                                \normalfont\bfseries #1}

\newenvironment{quotation}
               {\list{}{\listparindent 1.5em%
                        \itemindent    \listparindent
                        \rightmargin   \leftmargin
                        \parsep        \z@ \@plus\p@}%
                \item\relax}
               {\endlist}
\newenvironment{quote}
               {\list{}{\rightmargin\leftmargin}%
                \item\relax}
               {\endlist}

\setlength\arraycolsep{5\p@}
\setlength\tabcolsep{6\p@}
\renewcommand{\arraystretch}{1.1}
\setlength\arrayrulewidth{.4\p@}
\setlength\doublerulesep{2\p@}
\setlength\tabbingsep{\labelsep}
\skip\@mpfootins = \skip\footins
\setlength\fboxsep{3\p@}
\setlength\fboxrule{.4\p@}
\renewcommand \theequation {\@arabic\c@equation}
\newcounter{figure}
\renewcommand \thefigure {\@arabic\c@figure}
\def\fps@figure{tbp}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\fnum@figure{\figurename\nobreakspace\thefigure}
\newenvironment{figure}
               {\@float{figure}}
               {\end@float}
\newenvironment{figure*}
               {\@dblfloat{figure}}
               {\end@dblfloat}
\newcounter{table}
\renewcommand\thetable{\@arabic\c@table}
\def\fps@table{tbp}
\def\ftype@table{2}
\def\ext@table{lot}
\def\fnum@table{\tablename\nobreakspace\thetable}
\newenvironment{table}
               {\@float{table}}
               {\end@float}
\newenvironment{table*}
               {\@dblfloat{table}}
               {\end@dblfloat}
\newlength\abovecaptionskip
\newlength\belowcaptionskip
%\setlength\abovecaptionskip{8\p@}
\renewcommand{\abovecaptionskip}{0pt}
%\setlength\belowcaptionskip{0\p@}
\renewcommand{\belowcaptionskip}{4pt}
\long\def\@makecaption#1#2{%
  \vskip\abovecaptionskip
  \sbox\@tempboxa{#1: #2}%
  \ifdim \wd\@tempboxa >\hsize
    #1: #2\par
  \else
    \global \@minipagefalse
    \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
  \fi
  \vskip\belowcaptionskip}
\DeclareOldFontCommand{\rm}{\normalfont\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\normalfont\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\normalfont\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\normalfont\bfseries}{\mathbf}
\DeclareOldFontCommand{\it}{\normalfont\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\normalfont\slshape}{\@nomath\sl}
\DeclareOldFontCommand{\sc}{\normalfont\scshape}{\@nomath\sc}
\DeclareRobustCommand*\cal{\@fontswitch\relax\mathcal}
\DeclareRobustCommand*\mit{\@fontswitch\relax\mathnormal}
\newcommand\@pnumwidth{1.55em}
\newcommand\@tocrmarg{2.55em}
\newcommand\@dotsep{4.5}
\setcounter{tocdepth}{3}
\newcommand\tableofcontents{%
    \section*{\contentsname
        \@mkboth{%
           \MakeUppercase\contentsname}{\MakeUppercase\contentsname}}%
    \@starttoc{toc}%
    }
\newcommand*\l@part[2]{%
  \ifnum \c@tocdepth >-2\relax
    \addpenalty\@secpenalty
    \addvspace{2.25em \@plus\p@}%
    \setlength\@tempdima{3em}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      {\leavevmode
       \large \bfseries #1\hfil \hb@xt@\@pnumwidth{\hss #2}}\par
       \nobreak
       \if@compatibility
         \global\@nobreaktrue
         \everypar{\global\@nobreakfalse\everypar{}}%
      \fi
    \endgroup
  \fi}
\newcommand*\l@section[2]{%
  \ifnum \c@tocdepth >\z@
    \addpenalty\@secpenalty
    \addvspace{1.0em \@plus\p@}%
    \setlength\@tempdima{1.5em}%
    \begingroup
      \parindent \z@ \rightskip \@pnumwidth
      \parfillskip -\@pnumwidth
      \leavevmode \bfseries
      \advance\leftskip\@tempdima
      \hskip -\leftskip
      #1\nobreak\hfil \nobreak\hb@xt@\@pnumwidth{\hss #2}\par
    \endgroup
  \fi}
\newcommand*\l@subsection{\@dottedtocline{2}{1.5em}{2.3em}}
\newcommand*\l@subsubsection{\@dottedtocline{3}{3.8em}{3.2em}}
\newcommand*\l@paragraph{\@dottedtocline{4}{7.0em}{4.1em}}
\newcommand*\l@subparagraph{\@dottedtocline{5}{10em}{5em}}
\newcommand\listoffigures{%
    \section*{\listfigurename}%
      \@mkboth{\MakeUppercase\listfigurename}%
              {\MakeUppercase\listfigurename}%
    \@starttoc{lof}%
    }
\newcommand*\l@figure{\@dottedtocline{1}{1.5em}{2.3em}}
\newcommand\listoftables{%
    \section*{\listtablename}%
      \@mkboth{%
          \MakeUppercase\listtablename}%
         {\MakeUppercase\listtablename}%
    \@starttoc{lot}%
    }
\let\l@table\l@figure
\newdimen\bibindent
\setlength\bibindent{1.5em}
\newenvironment{thebibliography}[1]
     {\section*{\refname}%
      \@mkboth{\MakeUppercase\refname}{\MakeUppercase\refname}%
      \list{\@biblabel{\@arabic\c@enumiv}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
            \leftmargin\labelwidth
            \advance\leftmargin\labelsep
            \@openbib@code
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            \renewcommand\theenumiv{\@arabic\c@enumiv}}%
      \sloppy
      \clubpenalty10000
      \@clubpenalty \clubpenalty
      \widowpenalty10000%
      \sfcode`\.\@m}
     {\def\@noitemerr
       {\@latex@warning{Empty `thebibliography' environment}}%
      \endlist}
\newcommand\newblock{\hskip .11em\@plus.33em\@minus.07em}
\let\@openbib@code\@empty
\newenvironment{theindex}
               {\if@twocolumn
                  \@restonecolfalse
                \else
                  \@restonecoltrue
                \fi
                \twocolumn[\section*{\indexname}]%
                \@mkboth{\MakeUppercase\indexname}%
                        {\MakeUppercase\indexname}%
                \thispagestyle{plain}\parindent\z@
                \parskip\z@ \@plus .3\p@\relax
                \columnseprule \z@
                \columnsep 35\p@
                \let\item\@idxitem}
               {\if@restonecol\onecolumn\else\clearpage\fi}
\newcommand\@idxitem{\par\hangindent 40\p@}
\newcommand\subitem{\@idxitem \hspace*{20\p@}}
\newcommand\subsubitem{\@idxitem \hspace*{30\p@}}
\newcommand\indexspace{\par \vskip 10\p@ \@plus5\p@ \@minus3\p@\relax}
\renewcommand\footnoterule{%
  \kern-3\p@
  \hrule\@width.4\columnwidth
  \kern2.6\p@}
\newcommand\@makefntext[1]{%
    \parindent 1em%
    \noindent
    \hb@xt@1.8em{\hss\@makefnmark}#1}
\def\today{\ifcase\month\or
  January\or February\or March\or April\or May\or June\or
  July\or August\or September\or October\or November\or December\fi
  \space\number\day, \number\year}
\setlength\columnsep{10\p@}
\setlength\columnseprule{0\p@}
\pagestyle{plain}
\pagenumbering{arabic}
\if@twoside
\else
  \raggedbottom
\fi
\if@twocolumn
  \twocolumn
  \sloppy
  \flushbottom
\else
  \onecolumn
\fi

\RequirePackage[small,sc,hang]{caption}             % utilizado para dar formato a los epígrafes
\RequirePackage[abbr]{harvard}                      % utilizado para determinar el estilo de la bibliografia

\def\harvardurl#1{*#1}

\clubpenalty10000
\widowpenalty10000

\endinput

%%
%% End of file `revcoles.cls'.

\labeldocument[%
firstpage = 1,%
volume = 31, number = 2,%
month = 6, year = 2008,%
monthreceived = 1, yearreceived = 2008,%
monthaccepted = 1, yearaccepted = 2008]%

ubicar el membrete déspues del \documentclass...

Para títulos de apendice se utiliza el comando \appendixsection

author
referees
edition
publication
symposium
project


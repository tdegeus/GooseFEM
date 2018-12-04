% ------------------------------------------------------------------------------
%  @LaTeX-class-file{
%     filename        = "goose-article.cls",
%     version         = "0.1.1",
%     date            = "1 April 2018",
%     codetable       = "ISO/ASCII",
%     keywords        = "LaTeX, goose-article",
%     supported       = "https://github.com/tdegeus/goose-article",
%     docstring       = "A customized article-class."
% ------------------------------------------------------------------------------

% ==============================================================================
% Basic settings
% ==============================================================================

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{goose-article}[2016/11/16 v0.3 customized article class]

% load default article class
\LoadClass[a4paper,fleqn,twoside]{article}

% ==============================================================================
% Options
% ==============================================================================

% bibliography [namecite]
% -----------------------

\newif\if@namecite
\let\if@namecite\iffalse
\DeclareOption{namecite}{\let\if@namecite\iftrue}

% Review: [narrow,doublespacing]
% ------------------------------

\newif\if@narrow
\let\if@narrow\iffalse
\DeclareOption{narrow}{\let\if@narrow\iftrue}

\newif\if@doublespacing
\let\if@doublespacing\iffalse
\DeclareOption{doublespacing}{\let\if@doublespacing\iftrue}

\newif\if@twocolumn
\let\if@twocolumn\iffalse
\DeclareOption{twocolumn}{\let\if@twocolumn\iftrue}

% No header
% ---------

\newif\if@empty
\let\if@empty\iffalse
\DeclareOption{empty}{\let\if@empty\iftrue}

% Font: [times,garamond,verdana]
% ------------------------------

\newif\if@font
\let\if@font\iffalse

\newif\if@times
\let\if@times\iffalse
\DeclareOption{times}{\let\if@times\iftrue\let\if@font\iftrue}

\newif\if@garamond
\let\if@garamond\iffalse
\DeclareOption{garamond}{\let\if@garamond\iftrue\let\if@font\iftrue}

\newif\if@verdana
\let\if@verdana\iffalse
\DeclareOption{verdana}{\let\if@verdana\iftrue\let\if@font\iftrue}

% Process options
% ---------------

\ProcessOptions\relax

% ==============================================================================
% Packages
% ==============================================================================

\RequirePackage[left=25mm,right=25mm,top=30mm,bottom=30mm]{geometry}
\if@narrow
  \newgeometry{left=40mm,right=40mm,top=30mm,bottom=30mm}
\fi

\RequirePackage{fix-cm}
\RequirePackage{graphicx}
\RequirePackage[labelformat=simple]{subcaption}
\RequirePackage{amssymb,amsfonts,bm,cancel}
\RequirePackage{enumerate}
\RequirePackage[fleqn]{amsmath}
\RequirePackage[font={small},labelfont=bf,labelsep=period]{caption}
\RequirePackage[bf]{titlesec}
\RequirePackage[dvipsnames]{xcolor}
\RequirePackage{authblk}
\RequirePackage{hyperref}
\RequirePackage{multicol}

\renewcommand\thesubfigure{(\alph{subfigure})}

% ==============================================================================
% Citation style
% ==============================================================================

\if@namecite
  \RequirePackage{natbib}
  \let\oldbibliography\bibliography
  \renewcommand{\bibliography}[1]{%
    \bibliographystyle{apalike}
    \setlength{\bibsep}{3pt plus 0.3ex}
    \def\bibfont{\scriptsize}
    \oldbibliography{#1}
  }
\else
  \RequirePackage[square,sort&compress,numbers,comma]{natbib}
  \let\oldbibliography\bibliography
  \renewcommand{\bibliography}[1]{%
    \bibliographystyle{unsrtnat}
    \setlength{\bibsep}{3pt plus 0.3ex}
    \def\bibfont{\scriptsize}
    \begin{multicols}{2}\raggedright\oldbibliography{#1}\end{multicols}
  }
\fi

% define DOI command
\newcommand*{\doi}[1]{\sloppy\href{http://dx.doi.org/#1}{doi:~#1}}
\newcommand*{\arxivid}[1]{\sloppy\href{http://arxiv.org/abs/#1}{arXiv:~#1}}

% ==============================================================================
% Font
% ==============================================================================

\if@font
  \RequirePackage{mathspec}
  \RequirePackage{fontspec}
\fi

\if@garamond
  \AtBeginDocument{%
    \setmathfont(Digits)[Numbers=OldStyle]{Garamond}%
    \setromanfont[%
      BoldFont={Garamond Bold},%
      ItalicFont={Garamond Italic},%
      Mapping=tex-text%
    ]{Garamond}%
  }%
\fi

\if@verdana
  \AtBeginDocument{%
    \defaultfontfeatures{Mapping=tex-text}
    \setmainfont[Mapping=tex-text]{Verdana}
    \setsansfont{Verdana}
    \renewcommand*{\familydefault}{\sfdefault}
    \setmathfont(Greek,Digits,Latin){Verdana}
    \setmathrm{Verdana}
  }%
\fi

\if@times
  \AtBeginDocument{%
    \setmathfont(Digits)[Numbers=OldStyle]{Times New Roman}%
    \setromanfont[%
      BoldFont={Times New Roman Bold},%
      ItalicFont={Times New Roman Italic},%
      Mapping=tex-text%
    ]{Times New Roman}%
  }%
\fi

% ==============================================================================
% Spacing
% ==============================================================================

\if@doublespacing
  \AtBeginDocument{%
    \usepackage{setspace}
    \doublespacing
  }%
\fi

% ==============================================================================
% Two-columns
% ==============================================================================

\if@twocolumn
  \AtBeginDocument{%
    \twocolumn
  }%
\fi

% ==============================================================================
% Hyperref
% ==============================================================================

\AtBeginDocument{%
  \hypersetup{%
    pdftitle=\@title,
    citecolor=Blue,
    filecolor=black,
    linkcolor=Blue,
    urlcolor=Blue,
    breaklinks,
    hidelinks,
    bookmarksopen=true,
  }%
}

% ==============================================================================
% Headers
% ==============================================================================

% {command}{format}{label}{sep}{before-code}{after-code}
% - section
\titleformat
  {\section}
  {\normalfont\large\bfseries}
  {\thesection}
  {1em}
  {}
  {}
% - subsection
\titleformat
  {\subsection}
  {\normalfont\normalsize\bfseries}
  {\thesubsection}
  {1em}
  {}
  {}
% - paragraph
\titleformat
  {\paragraph}
  {\normalfont\normalsize\bfseries}
  {}
  {}
  {}
  {}

% spacing {left}{top}{bottom}, of:
% - section
\titlespacing{\section}
  {0pt}{1eM}{0.8em}
% - subsection
\titlespacing{\subsection}
  {0pt}{0.8em}{0.4em}
% - paragraph
\titlespacing{\paragraph}
  {0pt}{0.8em}{0em}

% no paragraph indent after header
\AtBeginDocument{%
  \makeatletter
  \let\orig@afterheading\@afterheading
  \def\@afterheading{%
     \@afterindentfalse
    \orig@afterheading}
  \makeatother
}

% ==============================================================================
% Title
% ==============================================================================

% title, authors, affiliation layout
% - contact
\newcommand{\contact}[1]{\date{{#1}}}
% - basic settings
\newcommand{\nl}{\vspace*{-.3em}\protect\\ \fontsize{9}{10.8}\itshape}
\renewcommand\Authfont{\fontsize{11}{14.4}\selectfont}
\renewcommand\Affilfont{\fontsize{9}{10.8}\itshape}
\renewcommand*{\Authsep}{, }
\renewcommand*{\Authand}{, }
\renewcommand*{\Authands}{, }
% - actual layout
\renewcommand{\maketitle}{%
  \newpage
  \null
  \vskip 2em%
  \begin{center}%
  \let \footnote \thanks
    {\Large\bfseries \@title \par}%
    \vskip 1.5em%
    {\normalsize
      \lineskip .5em%
      \begin{tabular}[t]{c}%
        \@author
      \end{tabular}\par}%
    \vskip 0.5em%
    {\small \@date}%
  \end{center}%
  \par
  \vskip 1.5em
  \thispagestyle{fancy}
}

% ==============================================================================
% Header and footer
% ==============================================================================

\RequirePackage{fancyhdr}
\setlength{\headheight}{16pt}
\pagestyle{fancy}
\fancyhead{}
\fancyfoot{}
\renewcommand{\headrulewidth}{0pt}
\renewcommand{\footrulewidth}{0pt}
\if@empty
\else
\fancyhead[RO,LE]{\thepage}
\fi
\newcommand{\header}[1]{\fancyhead[LO,RE]{#1}}
\newcommand{\headerOdd}[1]{\fancyhead[LO]{#1}}
\newcommand{\headerEven}[1]{\fancyhead[RE]{#1}}

% ==============================================================================
% Abstract and keywords
% ==============================================================================

\renewenvironment{abstract}{\section*{Abstract}\noindent\ignorespaces}{\ignorespacesafterend}
\newcommand{\keywords}[1]{\vspace*{0.5eM}\noindent\textbf{Keywords: }#1}

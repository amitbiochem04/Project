
####following thing need to be typed 
.Open Terminal and type nano .bash_profile
.Paste in the following lines:

export PS1="\[\033[36m\]\u\[\033[m\]@\[\033[32m\]\h:\[\033[33;1m\]\w\[\033[m\]\$ "
export CLICOLOR=1
export LSCOLORS=ExFxBxDxCxegedabagacad
alias ls='ls -GFh'
###then exit from the click follwoing command 

Hit Control+O to save, then Control+X to exit out of nano
#####oether way around 

export CLICOLOR=1
export LSCOLORS=GxFxCxDxBxegedabagaced
CLICOLOR=1 simply enables coloring of your terminal.

LSCOLORS=... specifies how to color specific items.

    Contact GitHub API Training Shop Blog About 



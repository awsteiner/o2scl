while ($_ = <STDIN>) {
        $_ =~ s/right: 0px;/left: 45em;/g;
        print $_
}



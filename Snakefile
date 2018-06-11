rule all:
    input:
        '{name}.hamburger'

rule make_hamburger:
    input:
        warm_bun='{name}.warm_bun',
        pickle='{name}.pickle',
        roasted_onion='{name}.roasted_onion',
        cheese='{name}.cheese',
        tomato='{name}.tomato',
        patty='{name}.patty',
        lettuce='{name}.lettuce'
    output:
        '{name}.hamburger'
    script:
        'hamburger/make/wrapper.py'

rule buy_tomato:
    output:
        temp('{name}.dirty_tomato')
    script:
        'tomato/buy/wrapper.py'

rule wash_tomato:
    input:
        dirty_tomato='{name}.dirty_tomato'
    output:
        temp('{name}.tomato')
    script:
        'tomato/wash/wrapper.py'

rule buy_cucumber:
    output:
        temp('{name}.cucumber')
    script:
        'cucumber/buy/wrapper.py'

rule buy_vinegar:
    output:
        temp('{name}.vinegar')
    script:
        'vinegar/buy/wrapper.py'

rule make_vinegar:
    input:
        cucumber='{name}.cucumber',
        vinegar='{name}.vinegar'
    output:
        temp('{name}.pickle')
    script:
        'pickle/make/wrapper.py'

rule buy_onion:
    output:
        temp('{name}.dirty_onion')
    script:
        'onion/buy/wrapper.py'

rule wash_onion:
    input:
        dirty_onion='{name}.dirty_onion'
    output:
        temp('{name}.onion')
    script:
        'onion/wash/wrapper.py'

rule roast_onion:
    input:
        onion='{name}.onion'
    output:
        temp('{name}.roasted_onion')
    script:
        'onion/roast/wrapper.py'
    
rule buy_cheese:
    output:
        temp('{name}.cheese')
    script:
        'cheese/buy/wrapper.py'

rule buy_bun:
    output:
        temp('{name}.cold_bun')
    script:
        'bun/buy/wrapper.py'

rule warm_bun:
    input:
        cold_bun='{name}.cold_bun'
    output:
        temp('{name}.warm_bun')
    script:
        'bun/warm/wrapper.py'

rule buy_lettuce:
    output:
        temp('{name}.dirty_lettuce')
    script:
        'lettuce/buy/wrapper.py'

rule wash_lettuce:
    input:
        dirty_lettuce='{name}.dirty_lettuce'
    output:
        temp('{name}.lettuce')
    script:
        'lettuce/wash/wrapper.py'

rule buy_beef:
    output:
        temp('{name}.beef')
    script:
        'beef/buy/wrapper.py'

rule make_patty:
    input:
        beef='{name}.beef'
    output:
        temp('{name}.patty')
    script:
        'patty/make/wrapper.py'

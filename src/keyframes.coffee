class @keyframes

    @frameCollection = []
    
    @cache = [0]
    
    @expando = 'data' + +new Date()
    
    @data = (elem) ->
        cacheIndex = elem[@expando]
        nextCacheIndex = @cache.length
 
        if(!cacheIndex)
            cacheIndex = elem[@expando] = nextCacheIndex
            @cache[cacheIndex] = {}
 
        return {
            get: (key) -> 
                return keyframes.cache[cacheIndex][key]
            set: (key, val) -> 
                keyframes.cache[cacheIndex][key] = val
                return val
        }
    
    @styleElem = ''
    
    @browserCode = ''
    
    @setBrowserCode = ->
        ua = navigator.userAgent
        if ua.indexOf('Opera') != -1
            return '-o-'
        else if ua.indexOf('MSIE') != -1
            return '-ms-'
        else if ua.indexOf('WebKit') != -1
            return '-webkit-'
        else if navigator.product == 'Gecko'
            return '-moz-'
        else
            return ''
        
    @generate = ->
        keyframes.styleElem.innerHTML = ''
        browserType = @browserCode
        list = []
        
        for key,frameName of @frameCollection
            css = '@'+browserType+'keyframes '+frameName.data.name+'{'
            
            for key,frameData of frameName.data
                if key != 'name'
                    css += key + '{'
                    css += frameData + ';}'
    
            css += '}\n'
            list.push(css)
    
        list.push(' .boostKeyframe{transform:scale3d(1,1,1);}')
    
        @styleElem.innerHTML = list.join('')
    
    @removeHead = ->
        @styleElem.innerHTML = ''

    @reset = (elem, callback) ->
        animationkey = @browserCode + 'animation'
        elem.style[animationkey + '-play-state'] = 'running'
        elem.style[animationkey] = 'none'
        @data(elem).set("keyframe", false)
        clearInterval(@data(elem).get('keyframeTimer'))
        clearTimeout(@data(elem).get('keyframeTimer'))
        if(callback)
            setTimeout(callback,0)
    
    @pause = (elem) ->
        elem.style[@browserCode + 'animation-play-state'] = 'paused'
        clearInterval(@data(elem).get('keyframeTimer'))
        clearTimeout(@data(elem).get('keyframeTimer'))
    
    @resume = (elem) ->
        elem.style[@browserCode + 'animation-play-state'] = 'running'
    
    @add = (frameData) ->
        for data in frameData
            kfname = data.name
            @frameCollection[kfname] = {data: data}
    
    @playNextItem = (elem, callback) ->
        frames = keyframes.data(elem).get('keyframeQueue')
        if frames.length
            keyframes.play(elem, frames[0], ->
                keyframes.playNextItem(elem, callback)
            )
            keyframes.data(elem).set('keyframeQueue', frames.slice(1))
        else
            if callback
                callback(elem)
    
    @play = (elem, frameOptions, callback, reset = true, regen = true) ->
        _playKeyframe = ->
            if typeof frameOptions == 'string'
                frameOptions = frameOptions.trim()
                frameOptSplit = frameOptions.split(' ')
                name = frameOptSplit[0]
                duration = parseInt(frameOptSplit[1])
                delay = parseInt(frameOptSplit[3])
                repeat = parseInt(frameOptSplit[4])
                frameOptSplit[1] += 'ms'
                frameOptSplit[3] += 'ms'
                
                animationcss = frameOptSplit.join(' ')
            else if frameOptions.length
                keyframes.data(elem).set('keyframeQueue', frameOptions)
                keyframes.playNextItem(elem, callback)
                return
            else
                name = frameOptions.name
                duration = frameOptions.duration
                delay = frameOptions.delay
                repeat = frameOptions.repeat
                frameOptions.duration += 'ms'
                frameOptions.delay += 'ms'
                animationcss = ''
                for key,opt of frameOptions
                    animationcss += opt + ' '
                    
                animationcss = animationcss.trim()
            
            animationkey = keyframes.browserCode + 'animation'
            
            if regen
                keyframes.generate()
            
            if repeat != 'infinite'
                if callback
                    keyframes.data(elem).set('keyframeTimer', setTimeout(->
                        callback(elem)
                    , (duration + delay) * repeat))
                setTimeout( ->
                    keyframes.data(elem).set('keyframe',false)
                , (duration + delay) * repeat)
            else
                if callback
                    keyframes.data(elem).set('keyframeTimer', setTimeout( ->
                        callback(elem)
                        keyframes.data(elem).set('keyframeTimer', setInterval(->
                            callback(elem)
                        , duration))
                    , duration + delay))
                
            
            elem.setAttribute("style", keyframes.browserCode + 'animation-play-state: running; '+animationkey + ': ' + animationcss)
            keyframes.data(elem).set('keyframe', name)
        
        if reset
            @reset(elem,_playKeyframe)
        else
            _playKeyframe()

document.getElementsByTagName('head')[0].innerHTML += '<style id="keyframes-style" type="text/css"></style>'
keyframes.styleElem = document.getElementById('keyframes-style')
keyframes.browserCode = keyframes.setBrowserCode()
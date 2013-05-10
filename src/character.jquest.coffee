class window.Character
    @settings: fileExtension: 'json', spriteDirectory: ''
    
    @sprites: []
    
    @characters: []
    
    @playerSprite: false
    
    @playerMoving: false
    
    @loadSprite: (charSource, callback) ->
        spritefile = "#{@settings.spriteDirectory+charSource}.#{@settings.fileExtension}"
        $.getJSON spritefile, (json) ->
            window.Character.sprites[json.id] = json
            if callback
                callback(json.id)
        
    @placeCharacter: (id, sprite, location) ->
        loc = window.Map.tileIdPosition(location)
        $("#layer#{Map.playerLayer}").append("<div class=\"character\" id=\"character#{id}\" />")
        charid = "#character#{id}"
        $(charid).width(@sprites[sprite].cellWidth).height(@sprites[sprite].cellHeight).css("background-image","url(#{@sprites[sprite].sprite})").css("background-position","-#{@sprites[sprite].cellWidth}px 0px").css("left",loc[0]).css("top",loc[1])
    
    @loadPlayer: (sprite) ->
        window.Character.playerSprite = sprite
        
        # ID 0 is always the player character.
        window.Character.characters[0] = sprite: sprite, timer: []
        window.Character.placeCharacter(0, sprite, window.Map.playerTile)
        
    @playerMove: (paths) ->
        ###
            If the character is already moving at the end location is redefined, it's best to remove the first
            of the next path as animation becomes more fluid.
        ###
        if @playerMoving
            paths = paths.slice(1)
            @stopAnimation(0)
            for index in [0...@characters[0].timer.length]
                clearTimeout(@characters[0].timer[index])
                
        $('#character0').clearQueue().stop()
        
        @playerMoving = true
        
        if paths[0]
            @animate(0,paths[0].direction)
            
        $.each paths, (index, path) ->
            loc = Map.tileIdPosition([path.x,path.y])
            
            if path.direction
                window.Character.characters[0].timer[index] = setTimeout -> 
                    window.Character.animate(0,path.direction)
                , index*310
            
            $('#character0').animate {
                left: loc[0]
                top: loc[1]
            },300,'linear', ->
                if index == paths.length-1
                    window.Character.playerMoving = false
                    window.Character.stopAnimation(0)
                Map.playerTile = Map.tileIdConvert([path.x+1,path.y])
            
    @animate: (id,direction) ->
        window.Character.stopAnimation(id)
        window.Character.characters[id].animation = {
            timer: setInterval ->
                if window.Character.characters[id].animation.step > 5
                    window.Character.characters[id].animation.step = 0
                y = 0
                switch (direction)
                    when "n"
                        y = 3
                    when "e"
                        y = 1
                    when "s"
                        y = 0
                    when "w"
                        y = 2
                    when "ne"
                        y = 6
                    when "nw"
                        y = 7
                    when "se"
                        y = 5
                    when "sw"
                        y = 4
                    
                y *= window.Character.sprites[id].cellHeight
                x = window.Character.characters[id].animation.step * window.Character.sprites[id].cellWidth
                $("#character#{id}").css("background-position", "-#{x}px -#{y}px")
                window.Character.characters[id].animation.step++
            , 200
        step: 0, y: 0 }
    
    @stopAnimation: (id) ->
        if window.Character.characters[id].animation
            clearInterval(window.Character.characters[id].animation.timer)
            bgpos = $("#character#{id}").css('background-position').split(' ')
            bgpos[0] = "-#{@sprites[id].cellWidth}px"
            $("#character#{id}").css('background-position', bgpos.join(' '))

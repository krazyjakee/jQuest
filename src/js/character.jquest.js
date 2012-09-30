var Character = {
	
	settings: {
		fileExtension: 'json',
		spriteDirectory: ''
	},
	sprites: [],
	characters: [],
	playerSprite: false,
	playerMoving: false,
	loadSprite: function(charSource, callback){
		var spritefile = Character.settings.spriteDirectory+charSource+'.'+Character.settings.fileExtension;
		$.getJSON(spritefile, function(json){
			Character.sprites[json.id] = json;
			if(callback){
				callback(json.id);
			}
		});
	},
	placeCharacter: function(id, sprite, location){
		var loc = Map.tileIdPosition(location);
		$('#layer'+Map.playerLayer).append('<div class="character" id="character'+id+'"></div>');
		
		var charid = '#character'+id;
		$(charid).width(Character.sprites[sprite].cellWidth).height(Character.sprites[sprite].cellHeight).css('background-image','url('+Character.sprites[sprite].sprite+')').css('background-position','-'+Character.sprites[sprite].cellWidth+'px 0px').css('left',loc[0]).css('top',loc[1]);
	},
	loadPlayer: function(sprite){
		Character.playerSprite = sprite;
		
		// ID 0 is always the player character.
		Character.characters[0] = { sprite: sprite, timer: [] };
		Character.placeCharacter(0, sprite, Map.playerTile);
	},
	playerMove: function(paths){
		
		/*
			If the character is already moving at the end location is redefined, it's best to remove the first
			of the next path as animation becomes more fluid.
		*/
		if(Character.playerMoving){
			paths = paths.slice(1);
			Character.stopAnimation(0);
			for(var i = 0; i <= Character.characters[0].timer.length; i++){
				clearTimeout(Character.characters[0].timer[i]);
			}
		}
		$('#character0').clearQueue().stop();
		
		Character.playerMoving = true;
		if(paths[0]){
			Character.animate(0,paths[0].direction);
		}
		$(paths).each(function(index,path){
			var loc = Map.tileIdPosition([path.x,path.y]);
			if(path.direction){
				Character.characters[0].timer[index] = setTimeout(function(){ Character.animate(0,path.direction); },index*310);
			}
			$('#character0').animate({
				left: loc[0],
				top: loc[1]
			},300,'linear',function(){
				if(index == paths.length-1){
					Character.playerMoving = false;
					Character.stopAnimation(0);
				}
				Map.playerTile = Map.tileIdConvert([path.x+1,path.y]);
			});
		});
		
	},
	animate: function(id,direction){
		Character.stopAnimation(id);
		
		Character.characters[id].animation = {timer: setInterval(function(){
			if(Character.characters[id].animation.step > 5){
				Character.characters[id].animation.step = 0;
			}
			var y = 0;
			switch(direction){
				case "n":
					var y = 3 * Character.sprites[id].cellHeight;
					break;
				case "e":
					var y = 1 * Character.sprites[id].cellHeight;
					break;
				case "s":
					var y = 0 * Character.sprites[id].cellHeight;
					break;
				case "w":
					var y = 2 * Character.sprites[id].cellHeight;
					break;
				case "ne":
					var y = 6 * Character.sprites[id].cellHeight;
					break;
				case "nw":
					var y = 7 * Character.sprites[id].cellHeight;
					break;
				case "se":
					var y = 5 * Character.sprites[id].cellHeight;
					break;
				case "sw":
					var y = 4 * Character.sprites[id].cellHeight;
					break;				
			}
			var x = Character.characters[id].animation.step * Character.sprites[id].cellWidth;
			$('#character'+id).css('background-position','-'+x+'px -'+y+'px');
			Character.characters[id].animation.step++
		},200), step: 0, y: 0 }
	},
	stopAnimation: function(id){
		if(Character.characters[id].animation){
			clearInterval(Character.characters[id].animation.timer);
			var bgpos = $('#character'+id).css('background-position').split(' ');
			bgpos[0] = '-'+Character.sprites[id].cellWidth+'px';
			$('#character'+id).css('background-position',bgpos.join(' '));
		}
	}
	
}
var Map = {
	settings: {
		fileExtension: 'json'
	},
	mapData: false,
	renderedTiles: [],
	tileProperties: [],
	playerTile: false,
	showBlocked: false,
	showPaths: false,
	playerLayer: 0,
	isMoving: false,
	drawMap: function(targetElem){
		$.each(Map.mapData.layers, function(index, layer){
			if(layer.type == "tilelayer"){
				$(targetElem).append('<div class="layer" id="layer'+index+'"></div>');
				
				if(layer.name == 'player'){
					Map.playerLayer = index;
				}
					
				for(var i = 0; i < layer.data.length; i++){
					$('#layer'+index).append('<div id="tile'+index+'-'+i+'" class="tile"></div>');
					Map.drawTile(index, layer.data[i], i);
				}
				
				$('#layer'+index).css('width',Map.mapData.tilewidth * Map.mapData.width + 'px');
				$('#layer'+index).css('height',Map.mapData.tileheight * Map.mapData.height + 'px');
				$('.tile').width(Map.mapData.tilewidth);
				$('.tile').height(Map.mapData.tileheight);
				$(targetElem).addClass('viewport');
				
				$('.layer:last').find('.tile').click(Map.tileClick);

				$('#maploading').remove();
			}else if(layer.type == "objectlayer"){
				
			}
		});
	},
	loadMap: function(targetElem, mapSource, callback){
		if(typeof mapSource == 'object'){
			Map.mapData = mapSource;
			Map.drawMap(targetElem);
			if(callback){
				callback();
			}
		}else if(typeof mapSource == 'string'){
			$.getJSON(mapSource+'.'+Map.settings.fileExtension, function(json){
				Map.mapData = json;
				Map.drawMap(targetElem);
				if(callback){
					callback();
				}
			});
		}
	},
	drawTile: function(layer, srcTile, targetTile){
		var target = $('#tile'+layer+'-'+targetTile);
			
		var setData = false;
		$.each(Map.mapData.tilesets, function(index, tileset){
			if(srcTile >= tileset.firstgid){
				setData = tileset;
			}else{
				return false;
			}
		});
		
		var setWidth = Math.floor(setData.imagewidth / setData.tilewidth)-1;
		var setHeight = Math.floor(setData.imageheight / setData.tileheight)-1;
		
		var heightCount = 0;
		var widthCount = 0;
		for(var i = 1; i < srcTile; i++){
			if(widthCount < setWidth){
				widthCount++;
			}else{
				widthCount = 0;
				heightCount++;
			}
		}
		
		var offset = {
			x: (widthCount * setData.tilewidth) + (setData.spacing * widthCount) + setData.margin,
			y: (heightCount * setData.tileheight) + (setData.spacing * heightCount) + setData.margin
		}
		
		if(setData.image){
			Map.renderedTiles[srcTile] = '#tile'+layer+'-'+targetTile;
			var property = Map.tileProperty(srcTile);
			if(property != false){
				Map.tileProperties[targetTile] = property;
			}
			$(target).css('background-image','url('+setData.image+')');
			$(target).css('background-position','-'+offset.x+'px -'+offset.y+'px');
		}
	},
	setFocus: function(tileId, duration){
		var loc = Map.tileIdConvert(tileId);
		var curloc = $('.layer:first').position();
		var maxwidth = (Map.mapData.width*32) - $('.viewport').width();
		var maxheight = (Map.mapData.height*32) - $('.viewport').height();
		loc[0] = (loc[0]*32) - ($('.viewport').width()/2);
		loc[1] = (loc[1]*32) - ($('.viewport').height()/2);
		
		// Set view limits
		if(loc[1] < 0){loc[1]=0;}
		if(loc[0] < 0){loc[0]=0;}
		if(loc[1] > maxheight){loc[1]=maxheight;}
		if(loc[0] > maxwidth){loc[0]=maxwidth;}
		
		Map.isMoving = true;
		
		$('.layer').clearQueue();
		$('.layer').animate({
	        top: "-"+loc[1]+"px",
			left: "-"+loc[0]+"px"
	    }, duration, function(){
	    	Map.isMoving = false;
	    });
	},
	tileIdConvert: function(tileInput){
		if(typeof tileInput == 'object'){
			var tileId = tileInput[1] * Map.mapData.width-1;
			tileId += tileInput[0];
			return tileId;
		}else{
			var y = 0;
			var x = 0;
			var width = Map.mapData.width-1;
			
			for(var i = 0; i < tileInput; i++){
				if(x < width){
					x++;
				}else{
					x = 0;
					y++;
				}
			}
			return [x, y];
		}
	},
	tileProperty: function(tileId){
		tileId = tileId - 1;
		var property = false;
		$.each(Map.mapData.tilesets, function(index, data){
			if(data.tileproperties != null && data.tileproperties[tileId] != null){
				property = data.tileproperties[tileId].property;
			}
		});
		return property;
	},
	tileClick: function(e){
		var tileId = $(this).attr('id');
		tileId = tileId.substr(tileId.lastIndexOf('-')+1);
		tileId++;
		if(e.button == 2){
			alert()
		}
		else{
			var paths = Map.makePath(tileId);
			if(paths.length){
				if(Map.showPaths){
					$('#layer'+Map.playerLayer+' div').css('background-color','transparent');
					$.each(paths, function(index, path){
						var tileId = Map.tileIdConvert([path.x,path.y]);
						$('#tile'+Map.playerLayer+'-'+tileId).css('background-color','red');
					});
				}
				Map.setFocus(tileId, (paths.length * 500));
				Map.playerTile = tileId;
			}
		}
	},
	makePath: function(toTileId){
		var totalMapSize = Map.mapData.width * Map.mapData.height;
		var toTileLoc = Map.tileIdConvert(toTileId);
		var fromTileLoc = Map.tileIdConvert(Map.playerTile);
		var board = [];
		
		for(var y = 0; y < Map.mapData.height-1; y++)
		{
			board[y] = [];
			for (var x = 0; x < Map.mapData.width-1; x++){
				var tile = Map.tileIdConvert([x,y]);
				var prop = Map.tileProperties[tile];
				if(prop == 'block'){
					board[y][x] = 1;
					if(Map.showBlocked){
						var tileidc = this.tileIdConvert([x,y]);
						$('#tile0-'+tileidc).css('background','red');
					}
				}else{
					board[y][x] = 0;
				}
			}
		}
		return AStar(board, fromTileLoc, toTileLoc, 'Diagonal');
	}
}
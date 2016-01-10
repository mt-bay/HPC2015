//------------------------------------------------------------------------------
/// @file
/// @brief    HPCAnswer.hpp の実装 (解答記述用ファイル)
/// @author   ハル研究所プログラミングコンテスト実行委員会
///
/// @copyright  Copyright (c) 2015 HAL Laboratory, Inc.
/// @attention  このファイルの利用は、同梱のREADMEにある
///             利用条件に従ってください

//------------------------------------------------------------------------------

#include "HPCAnswer.hpp"
#include "HPCMath.hpp"

#include <vector>
#include <list>
#include <map>
#include <algorithm>
#include <iterator>
#include <limits>
#include <cstdint>

//#include <iostream>

/// プロコン問題環境を表します。
namespace hpc {
    int distance(const Pos& _a, const Pos& _b);

    void choose_item(const Field& _field, const ItemCollection& _items);
    std::vector<Action> gen_route(const Field& _field, const ItemCollection& _items, const std::list<int>& _load_index);
    std::list<int> convex_hull(const Field& _field, const ItemCollection& _items, const std::vector<int>& _load_index);
    std::vector<int> sort_position_y(const Field& _field, const ItemCollection& _items, const std::vector<int>& _target);
    const Pos get_pos(const Field& _field, const ItemCollection& _items, const int _index);
    int distance(const Pos& _a, const Pos& _b);

    //------------------------------------------------------------------------------
    /// @brief startからendまでの経路を示します
    class Path
    {
    public:
        int start;
        int goal;

        std::list<Action> route;

    private:
        class AStar_dat
        {
        private:
            inline int get_huries(const Pos& _goal, const Pos& _here)
                const
            {
                return Math::Abs(_goal.x - _here.x) + Math::Abs(_goal.y - _here.y);
            }

        public:
            enum State
            {
                NONE,
                OPEN,
                CLOSE
            };

            AStar_dat*    m_parent;
            Pos   m_pos;

            State m_state;
            int   m_cost;
            int   m_huries;

            // デフォルトコンストラクタ
            AStar_dat()
                :
                m_parent(nullptr),
                m_state (NONE)   ,
                m_cost  (0)      ,
                m_huries(0)
            {
            }

            // コピーコンストラクタ
            AStar_dat(const AStar_dat& _orig)
                :
                m_parent(_orig.m_parent),
                m_state (_orig.m_state) ,
                m_cost  (_orig.m_cost)  ,
                m_huries(_orig.m_huries)
            {
            }

            // ノードを開く
            // @return ノードがOPENになったか
            bool open(const Pos& _here, const Pos& _goal, AStar_dat* _parent = nullptr, const bool _is_wall = false)
            {
                if(_is_wall)
                {
                    m_state = CLOSE;
                    return false;
                }

                m_state  = OPEN;
                m_parent = _parent;
                m_pos    = _here;
                m_cost   = (m_parent == nullptr ? 0 : m_parent->m_cost) + 1;
                m_huries  = get_huries(_goal, _here);

                return true;
            }

            inline int get_score()
                const
            {
                return m_cost + m_huries;
            }
        };

        // ゴールからスタートまでのノードがOPENになっているデータからPATH導出
        std::list<Action> gen_path(std::vector< std::vector<AStar_dat> >& _data, const Pos& _goal)
            const
        {
            std::list<Action> buf_return;

            AStar_dat* buf_data = &_data.at(_goal.x).at(_goal.y), *buf_parent = buf_data->m_parent;
            while(buf_parent != nullptr)
            {
                if     (buf_data->m_pos.x - buf_parent->m_pos.x ==  1)
                {
                    buf_return.push_front(Action_MoveRight);
                }
                else if(buf_data->m_pos.x - buf_parent->m_pos.x == -1)
                {
                    buf_return.push_front(Action_MoveLeft);
                }
                else if(buf_data->m_pos.y - buf_parent->m_pos.y ==  1)
                {
                    buf_return.push_front(Action_MoveUp);
                }
                else if(buf_data->m_pos.y - buf_parent->m_pos.y == -1)
                {
                    buf_return.push_front(Action_MoveDown);
                }
                buf_data   = &_data.at(buf_parent->m_pos.x).at(buf_parent->m_pos.y);
                buf_parent = buf_parent->m_parent;
            }

            return buf_return;
        }

        // A*探索
        std::list<Action> a_star(const Field& _field, const Pos& _start, const Pos& _goal)
            const
        {
            if(_goal == _start)
            {
                return std::list<Action>();
            }

            std::vector< std::vector<AStar_dat> > data(_field.width(), std::vector<AStar_dat>(_field.height(), AStar_dat()));

            std::vector<AStar_dat*> open_nodes;  //Open状態のノードリスト
            std::list<AStar_dat*>   close_nodes; //Close状態のノードリスト

            //初期ノード生成
            data.at(_start.x).at(_start.y).open(_start, _goal);
            open_nodes.push_back(&data.at(_start.x).at(_start.y));

            Pos    buf_pos, buf_move;

            AStar_dat* ref_buf_data;
            //探索
            int open_num = 0;
            while(true)
            {
                // スコア最小のノードを導出
                int buf_index = 0, buf_score = open_nodes.front()->get_score();
                for(size_t s = 1, cond_s = open_nodes.size(); s < cond_s; ++s)
                {
                    if(buf_score > open_nodes.at(s)->get_score())
                    {
                        buf_index = s;
                        buf_score = open_nodes.at(s)->get_score();
                    }
                }
                buf_pos = open_nodes.at(buf_index)->m_pos;
                // スコア最小のノードの上下左右を探索する(OPENにする)
                for(int i = 0; i < 4; ++i)
                {
                    buf_move = buf_pos.move((Action)i);
                    ref_buf_data = &data.at(buf_move.x).at(buf_move.y);
                    if(ref_buf_data->m_state == AStar_dat::NONE)
                    {
                        if(ref_buf_data->open(buf_move, _goal, open_nodes.at(buf_index), _field.isWall(buf_move)))
                        {
                            open_nodes.push_back(ref_buf_data);
                        }
                        ++open_num;
                        if(ref_buf_data->m_pos == _goal)
                        {
                            return gen_path(data, _goal);
                        }
                    }
                }
                open_nodes.at(buf_index)->m_state = AStar_dat::CLOSE;
                open_nodes.erase(open_nodes.begin() + buf_index);
            }
        }

    public:
        // A*を用いてstartからgoalまでの経路導出
        Path(const Field& _field ,const ItemCollection& _item ,const int _start, const int _goal)
            :
            start(_start),
            goal (_goal),
            route(a_star(_field, get_pos(_field, _item, _start), get_pos(_field, _item, _goal)))
        {
        }

        Path()
            :
            start(0),
            goal (0),
            route(std::list<Action>())
        {
        }

        //コピーコンストラクタ
        Path(const Path& _orig)
            :
            start(_orig.start),
            goal (_orig.goal) ,
            route(_orig.route)
        {

        }

        Path& operator=(const Path& _rhs)
        {
            start = _rhs.start;
            goal  = _rhs.goal;
            route = _rhs.route;

            return *this;
        }

        inline size_t cost()
            const
        {
            return route.size();
        }

        inline std::list<Action> get_route()
            const
        {
            return route;
        }

        inline bool start_is_office(const int item_count)
            const
        {
            return start < 0 || item_count <= start;
        }
    };

    namespace 
    {
        std::vector<Action> action[Parameter::PeriodCount];
        std::list<int>      load_index[Parameter::PeriodCount];
        int    period;
        size_t act_index;
    };


    //------------------------------------------------------------------------------
    /// 各ステージ開始時に呼び出されます。
    ///
    /// ここで、各ステージに対して初期処理を行うことができます。
    ///
    /// @param[in] aStage 現在のステージ。
    void Answer::Init(const Stage& aStage)
    {
        choose_item(aStage.field(), aStage.items());

        for(size_t s = 0; s < Parameter::PeriodCount; ++s)
        {
            action[s] = gen_route(aStage.field(), aStage.items(), load_index[s]);
        }
    }

    //------------------------------------------------------------------------------
    /// 各配達時間帯開始時に呼び出されます。
    ///
    /// ここで、この時間帯に配達する荷物をトラックに積み込みます。
    /// どの荷物をトラックに積み込むかを決めて、引数の aItemGroup に対して
    /// setItem で荷物番号を指定して積み込んでください。
    ///
    /// @param[in] aStage 現在のステージ。
    /// @param[in] aItemGroup 荷物グループ。
    void Answer::InitPeriod(const Stage& aStage, ItemGroup& aItemGroup)
    {
        for(std::list<int>::const_iterator it = load_index[aStage.period()].cbegin(), cond_it = load_index[aStage.period()].cend(); it != cond_it; ++it)
        {
            aItemGroup.addItem(*it);
        }
        
        period    = aStage.period();
        act_index = 0;
        return;

        if (aStage.period() == 0) {
            return;
        }
        for (int i = 0; i < aStage.items().count(); ++i) {
            // まだ配達されてない荷物かどうか調べる
            if (aStage.getTransportState(i) == TransportState_NotTransported) {
                // 配達されてない荷物なので積み込む
                aItemGroup.addItem(i);
            }
        }
    }

    //------------------------------------------------------------------------------
    /// 各ターンでの動作を返します。
    ///
    /// @param[in] aStage 現在ステージの情報。
    ///
    /// @return これから行う動作を表す Action クラス。
    Action Answer::GetNextAction(const Stage& aStage)
    {
        return action[period].at(act_index++);

        static Random random; // デフォルトのシード値を使う
        static Pos prev; // 初期値は重要ではない。(前のゲームの値が残っていても気にしない)

        for (int retry = 0; ; ++retry) {
            Action a = static_cast<Action>(random.randTerm(4));
            Pos nextPos = aStage.truck().pos().move(a);
            if (aStage.field().isWall(nextPos) == false) { // 動けるか
                if (retry < 50 && nextPos == prev) {
                    // 前にいた場所を避ける。
                    // これで、同じような場所をウロウロしてなかなか進まないのを防げる。
                    // ただし、50回やっても見つからないときは、諦める。
                    continue;
                }
                prev = aStage.truck().pos();
                return a;
            }
        }
    }

    //------------------------------------------------------------------------------
    /// 各配達時間帯終了時に呼び出されます。
    ///
    /// ここで、この時間帯の終了処理を行うことができます。
    ///
    /// @param[in] aStage 現在のステージ。
    /// @param[in] aStageState 結果。Playingならこの時間帯の配達完了で、それ以外なら、何らかのエラーが発生した。
    /// @param[in] aCost この時間帯に消費した燃料。エラーなら0。
    void Answer::FinalizePeriod(const Stage& aStage, StageState aStageState, int aCost)
    {
        if (aStageState == StageState_Failed) {
            // 失敗したかどうかは、ここで検知できます。
        }
        else
        {
        }
    }

    //------------------------------------------------------------------------------
    /// 各ステージ終了時に呼び出されます。
    ///
    /// ここで、各ステージに対して終了処理を行うことができます。
    ///
    /// @param[in] aStage 現在のステージ。
    /// @param[in] aStageState 結果。Completeなら配達完了で、それ以外なら、何らかのエラーが発生した。
    /// @param[in] aScore このステージで獲得したスコア。エラーなら0。
    void Answer::Finalize(const Stage& aStage, StageState aStageState, int aScore)
    {
        if (aStageState == StageState_Failed) {
            // 失敗したかどうかは、ここで検知できます。
        }
        else if (aStageState == StageState_TurnLimit) {
            // ターン数オーバーしたかどうかは、ここで検知できます。
        }
        for(size_t s = 0; s < Parameter::PeriodCount; ++s)
        {
            load_index[s].clear();
            action[s].clear();
        }
    }


    //------------------------------------------------------------------------------
    //以下，追加関数
    
    //積載物選択
    //4ピリオドの荷物全部決定する
    void choose_item(const Field& _field, const ItemCollection& _items)
    {
        std::list<int> free_item[Parameter::ItemWeightMax];
        int period_weight[Parameter::PeriodCount] = {};
        //各ピリオドで配達することを指定されているものはとりあえず載せる
        for(size_t s = 0, cond_s = (size_t)_items.count(); s < cond_s; ++s)
        {
            if(_items[s].period() == -1)
            {
                free_item[_items[s].weight() - 1].push_back(s);
            }
            else
            {
                load_index[_items[s].period()].push_back(s);
                period_weight[_items[s].period()] += _items[s].weight();
            }
        }
        
        //時間帯指定なしの荷物を，重いものから載せていく
        //載せるピリオドが決まっているものから最近傍点を導出し，その点と同じピリオドに載せることにする
        int min_period, min_distance, buf_distance, target_weight = Parameter::ItemWeightMax;
        for(int i = Parameter::ItemWeightMax - 1; i >= 0; --i)
        {
            while(!free_item[i].empty())
            {
                min_period = 0;
                min_distance = std::numeric_limits<int32_t>::max();
                for(size_t t = 0; t < Parameter::PeriodCount; ++t)
                {
                    if(target_weight + period_weight[t] > Parameter::TruckWeightCapacity)
                    {
                        continue;
                    }

                    if(load_index[t].empty())
                    {
                        if(min_distance > (buf_distance = distance(get_pos(_field, _items, free_item[i].front()), _field.officePos())))
                        {
                            min_period = t;
                            min_distance = buf_distance;
                        }
                    }

                    for(std::list<int>::const_iterator it = load_index[t].cbegin();
                        it != load_index[t].cend();
                        ++it)
                    {
                        if(min_distance > (buf_distance = distance(get_pos(_field, _items, free_item[i].front()), get_pos(_field, _items, *it))))
                        {
                            min_period = t;
                            min_distance = buf_distance;
                        }
                    }
                }
                

                load_index[min_period].push_back(free_item[i].front());
                period_weight[min_period] += target_weight;
                free_item[i].pop_front();
            }
            --target_weight;
        }

        //営業所を表すindexを全てのピリオドに追加
        for(size_t s = 0; s < Parameter::PeriodCount; ++s)
        {
            load_index[s].unique();
        }
    }

    //巡回ルート生成
    std::vector<Action> gen_route(const Field& _field, const ItemCollection& _items, const std::list<int>& _load_index)
    {
        if(_load_index.empty())
        {
            return std::vector<Action>();
        }

        std::list<int> load_index(_load_index);
        load_index.push_back(_items.count());

        //目的地集合の凸包(エッジは壁を無視して考える)を作り目的地の巡回ルートとする
        std::list<Path> buf_route;
        std::list<int>  buf_dist(convex_hull(_field, _items, std::vector<int>(load_index.begin(), load_index.end())));
        std::list<int>  last;

        for(std::list<int>::const_iterator s = load_index.cbegin(), cond_s = load_index.cend(); s != cond_s; ++s)
        {
            if(std::find(buf_dist.begin(), buf_dist.end(), *s) == buf_dist.end())
            {
                last.push_back(*s);
            }
        }

        //巡回ルートの経路を導出する
        std::vector< std::map<size_t, Path> > cache_route(_items.count() + 1, std::map<size_t, Path>());
        std::vector< std::map<size_t, int>  > cache_cost (_items.count() + 1, std::map<size_t, int>() );
        for(std::list<int>::const_iterator it_a = buf_dist.cbegin(), it_b = ++buf_dist.cbegin(), cond_it = buf_dist.cend(); it_b != cond_it; )
        {
            buf_route.push_back(Path(_field, _items, *it_a, *it_b));
            cache_route.at(*it_a).insert(std::map<size_t, Path>::value_type(*it_b, buf_route.back()));
            cache_cost.at(*it_a).insert(std::map<size_t, int>::value_type(*it_b, buf_route.back().cost()));
            ++it_a; ++it_b;
        }
        buf_route.push_back(Path(_field, _items, buf_dist.back(), buf_dist.front()));
        cache_route.at(buf_dist.back()).insert(std::map<size_t, Path>::value_type(buf_dist.front(), buf_route.back()));
        cache_cost.at(buf_dist.back()).insert(std::map<size_t, int>::value_type(buf_dist.front(), buf_route.back().cost()));

        //元々のルートから目的地に選択されていないノードを経由させた場合の距離差分が最小となるものに貼り直す
        // これを全目的地が巡回ルートに入るように繰り返す
        std::list<int>::iterator min_goal, min_last;
        int min_diff;
        
        std::list<Path>::iterator it_route_a = buf_route.begin(), min_route_start;
        while(!last.empty())
        {
            it_route_a = buf_route.begin();
            min_diff = std::numeric_limits<std::int32_t>::max();
            for(std::list<int>::iterator it_dist_a = buf_dist.begin(), it_dist_b = ++buf_dist.begin(), cond_dist_it = buf_dist.end(); it_dist_b != cond_dist_it; )
            {
                for(std::list<int>::iterator it_last = last.begin(), cond_it_last = last.end(); it_last != cond_it_last; ++it_last)
                {
                    if(cache_cost.at(*it_dist_a).count(*it_last) == 0)
                    {
                        cache_route.at(*it_dist_a)[*it_last] = Path(_field, _items, *it_dist_a, *it_last);
                        cache_cost .at(*it_dist_a)[*it_last] = cache_route.at(*it_dist_a)[*it_last].cost();
                    }
                    if(cache_cost.at(*it_last).count(*it_dist_b) == 0)
                    {
                        cache_route.at(*it_last)[*it_dist_b] = Path(_field, _items, *it_last, *it_dist_b);
                        cache_cost .at(*it_last)[*it_dist_b] = cache_route.at(*it_last)[*it_dist_b].cost();
                    }
                    if(min_diff > cache_cost.at(*it_dist_a).at(*it_last) + cache_cost.at(*it_last).at(*it_dist_b) - cache_cost.at(*it_dist_a).at(*it_dist_b))
                    {
                        min_route_start = it_route_a;
                        min_goal = it_dist_b; min_last = it_last;
                        min_diff = cache_cost.at(*it_dist_a).at(*it_last) + cache_cost.at(*it_last).at(*it_dist_b) - cache_cost.at(*it_dist_a).at(*it_dist_b);
                    }
                }
                ++it_dist_a;
                ++it_dist_b;
                ++it_route_a;
            }

            std::list<int>::iterator min_cp    = buf_dist.insert(min_goal, *min_last);
            std::list<int>::iterator min_start = std::next(min_cp, -1);
            min_goal = std::next(min_cp, 1);

            min_route_start = buf_route.erase(min_route_start);
            min_route_start = buf_route.insert(min_route_start, cache_route.at(*min_cp).at(*min_goal));
            min_route_start = buf_route.insert(min_route_start, cache_route.at(*min_start).at(*min_cp));

            last.erase(min_last);
        }
        //始点を営業所とするように，アクションを結合
        std::list<Action> act_front, act_back;
        std::list<Path>::iterator it = buf_route.begin(), cond_it = buf_route.end();
        for(; it != cond_it; ++it)
        {
            if(it->start_is_office(_items.count()))
            {
                break;
            }
            act_back.insert(act_back.end(), it->route.begin(), it->route.end());
        }
        for(; it != cond_it; ++it)
        {
            act_front.insert(act_front.end(), it->route.begin(), it->route.end());
        }
        
        act_front.insert(act_front.end(), act_back.begin(), act_back.end());

        return std::vector<Action>(act_front.begin(), act_front.end());
    }


    std::list<int> convex_hull(const Field& _field, const ItemCollection& _items, const std::vector<int>& _load_index)
    {
        if(_load_index.size() <= 3)
        {
            return std::list<int>(_load_index.begin(), _load_index.end());
        }

        std::list<int> hull_index;
        //y軸最小のものを凸包の結果に挿入
        std::vector<int> sorted_load(sort_position_y(_field, _items, _load_index));

        hull_index.push_back(sorted_load.front());

        //角度から凸包を作る
        int min_index = 0;
        float min_angle = std::numeric_limits<float>::max(), buf_angle;
        //上から下
        for(int i = min_index, j, cond_i = sorted_load.size() - 1, cond_j = sorted_load.size(); i < cond_i; i = min_index)
        {
            min_angle = std::numeric_limits<float>::max();
            for(j = i + 1; j < cond_j; ++j)
            {
                if(min_angle > (buf_angle = Math::ATan2((float)(get_pos(_field, _items, sorted_load.at(i)).y - get_pos(_field, _items, sorted_load.at(j)).y),
                                                        (float)(get_pos(_field, _items, sorted_load.at(i)).x - get_pos(_field, _items, sorted_load.at(j)).x))))
                {
                    min_angle = buf_angle;
                    min_index = j;
                }
            }

            hull_index.push_back(sorted_load.at(min_index));
        }
        //下から上
        min_index = sorted_load.size() - 1;
        for(int i = min_index, j, cond_i = 0, cond_j = -1; i > cond_i; i = min_index)
        {
            min_angle = std::numeric_limits<float>::max();
            for(j = i - 1; j > cond_j; --j)
            {
                if(min_angle > (buf_angle = Math::ATan2((float)(get_pos(_field, _items, sorted_load.at(i)).y - get_pos(_field, _items, sorted_load.at(j)).y),
                                                        (float)(get_pos(_field, _items, sorted_load.at(i)).x - get_pos(_field, _items, sorted_load.at(j)).x))))
                {
                    min_angle = buf_angle;
                    min_index = j;
                }
            }
            hull_index.push_back(sorted_load.at(min_index));
        }

        //最初と最後の要素は被るハズなので最後の要素を削除
        hull_index.pop_back();

        return hull_index;
    }


    //y軸で参照indexを挿入ソート
    std::vector<int> sort_position_y(const Field& _field, const ItemCollection& _items, const std::vector<int>& _target)
    {
        if(_target.size() <= 1)
        {
            return _target;
        }
        std::vector<int> buf_return(_target);
        int buf;

        for(size_t s = 1, t, cond_s = buf_return.size(); s < cond_s; ++s)
        {
            for(t = 0; t < s; ++t)
            {
                if(get_pos(_field, _items, buf_return.at(t)).y > get_pos(_field, _items, buf_return.at(s)).y)
                {
                    break;
                }
            }
            for(; t < s; ++t)
            {
                buf = buf_return.at(t);
                buf_return.at(t) = buf_return.at(s);
                buf_return.at(s) = buf;
            }
        }

        return buf_return;
    }


    const Pos get_pos(const Field& _field, const ItemCollection& _items, const int _index)
    {
        return (int)0 <= _index && _index < _items.count() ? _items[_index].destination() : _field.officePos();
    }

    int distance(const Pos& _a, const Pos& _b)
    {
        return Math::Abs(_b.x - _a.x) + Math::Abs(_b.y - _a.y);
    }
}

//------------------------------------------------------------------------------
// EOF
